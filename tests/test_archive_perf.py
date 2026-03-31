"""
Tests for the optimised archive / merge_archive / unarchive methods.

Covers:
- archive round-trip (compress + uncompressed)
- merge_archive from compressed tarball (streaming extraction)
- merge_archive from uncompressed directory (parallel copy)
- CMSO.hasPath triples are rewritten correctly after each operation
- Structure-store files are present after each operation
- Duplicate file warning in merge_archive
"""

import os
import shutil
import tempfile
import tarfile
import warnings

import pytest

from atomrdf import KnowledgeGraph
from atomrdf.namespace import CMSO
import atomrdf.build as build


# ── helpers ──────────────────────────────────────────────────────────────────


@pytest.fixture()
def workdir():
    """Run each test inside a fresh temporary directory, then clean up."""
    d = tempfile.mkdtemp(prefix="atomrdf_archive_test_")
    prev = os.getcwd()
    os.chdir(d)
    yield d
    os.chdir(prev)
    shutil.rmtree(d, ignore_errors=True)


def _make_kg(n_structures=3):
    """Build an in-memory KG with *n_structures* bulk structures."""
    kg = KnowledgeGraph()
    elements = ["Fe", "Cu", "Al", "Ni", "Ti", "Cr", "Si", "Mg"]
    for i in range(n_structures):
        build.bulk(elements[i % len(elements)], graph=kg)
    return kg


def _has_path_basenames(kg):
    """Return set of basenames referenced by CMSO.hasPath triples."""
    basenames = set()
    for _, _, obj in kg.triples((None, CMSO.hasPath, None)):
        try:
            basenames.add(os.path.basename(obj.toPython()))
        except Exception:
            pass
    return basenames


def _structure_store_files(kg):
    """Return set of filenames present in the KG's structure store dir."""
    if os.path.isdir(kg.structure_store):
        return set(os.listdir(kg.structure_store))
    return set()


# ── archive + unarchive (compressed) ────────────────────────────────────────


class TestArchiveCompressed:

    def test_roundtrip_preserves_samples(self, workdir):
        kg = _make_kg(4)
        kg.archive("pkg_comp")
        assert os.path.exists("pkg_comp.tar.gz")

        kg2 = KnowledgeGraph.unarchive("pkg_comp.tar.gz")
        assert kg2.n_samples == 4

    def test_structure_files_inside_tarball(self, workdir):
        kg = _make_kg(2)
        kg.archive("pkg_files")

        with tarfile.open("pkg_files.tar.gz") as tar:
            names = tar.getnames()
        ss_files = [
            n for n in names if "rdf_structure_store" in n and not n.endswith("/")
        ]
        # At least one JSON per structure
        assert len(ss_files) >= 2

    def test_has_path_rewritten_in_archive(self, workdir):
        """After archive, hasPath triples should use rdf_structure_store/ prefix."""
        kg = _make_kg(2)
        kg.archive("pkg_paths", compress=False)

        kg2 = KnowledgeGraph(
            graph_file="pkg_paths/triples",
            structure_store="pkg_paths/rdf_structure_store",
        )
        for _, _, obj in kg2.triples((None, CMSO.hasPath, None)):
            path = obj.toPython()
            assert "rdf_structure_store/" in path

    def test_compresslevel_produces_valid_gzip(self, workdir):
        """The low-compresslevel tarball is still a valid gzip archive."""
        kg = _make_kg(2)
        kg.archive("pkg_gz")
        # tarfile.open will raise if the file is corrupt
        with tarfile.open("pkg_gz.tar.gz", "r:gz") as tar:
            assert len(tar.getnames()) > 0


# ── archive + unarchive (uncompressed) ──────────────────────────────────────


class TestArchiveUncompressed:

    def test_roundtrip_uncompressed(self, workdir):
        kg = _make_kg(3)
        kg.archive("pkg_unc", compress=False)
        assert os.path.isdir("pkg_unc")
        assert os.path.isdir("pkg_unc/rdf_structure_store")

        kg2 = KnowledgeGraph.unarchive("pkg_unc", compress=False)
        assert kg2.n_samples == 3

    def test_structure_store_files_copied(self, workdir):
        kg = _make_kg(2)
        kg.archive("pkg_unc2", compress=False)
        files = os.listdir("pkg_unc2/rdf_structure_store")
        assert len(files) >= 2


# ── merge_archive (compressed) ──────────────────────────────────────────────


class TestMergeArchiveCompressed:

    def test_merge_adds_samples(self, workdir):
        # Create two separate archives
        kg1 = _make_kg(2)
        kg1.archive("ds1")
        kg2 = _make_kg(3)
        kg2.archive("ds2")

        # Merge into a fresh KG
        merged = KnowledgeGraph()
        merged.merge_archive("ds1.tar.gz")
        assert merged.n_samples == 2
        merged.merge_archive("ds2.tar.gz")
        assert merged.n_samples == 5

    def test_streaming_extraction_copies_files(self, workdir):
        """Structure-store files should be streamed directly to structure_store."""
        kg = _make_kg(3)
        kg.archive("ds_stream")

        merged = KnowledgeGraph()
        merged.merge_archive("ds_stream.tar.gz")

        # Verify files exist in the merged KG's structure store
        ss_files = _structure_store_files(merged)
        path_basenames = _has_path_basenames(merged)
        # Every referenced path should have a corresponding file
        for bn in path_basenames:
            assert bn in ss_files, f"{bn} referenced in triple but missing from store"

    def test_has_path_rewritten_after_merge(self, workdir):
        kg = _make_kg(2)
        kg.archive("ds_rewrite")

        merged = KnowledgeGraph()
        merged.merge_archive("ds_rewrite.tar.gz")

        for _, _, obj in merged.triples((None, CMSO.hasPath, None)):
            path = obj.toPython()
            resolved = os.path.join(merged.structure_store, os.path.basename(path))
            assert os.path.exists(resolved), f"path {path} does not resolve to a file"

    def test_duplicate_file_warns(self, workdir):
        """Merging the same archive twice should warn about existing files."""
        kg = _make_kg(1)
        kg.archive("ds_dup")

        merged = KnowledgeGraph()
        merged.merge_archive("ds_dup.tar.gz")

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            merged.merge_archive("ds_dup.tar.gz")
            dup_warnings = [x for x in w if "already exists" in str(x.message)]
            assert len(dup_warnings) > 0

    def test_no_leftover_extracted_dir(self, workdir):
        """After merge_archive(compress=True), the extracted folder is cleaned up."""
        kg = _make_kg(2)
        kg.archive("ds_clean")

        merged = KnowledgeGraph()
        merged.merge_archive("ds_clean.tar.gz")

        # The extracted dir (ds_clean/) should have been removed
        assert not os.path.isdir("ds_clean")


# ── merge_archive (uncompressed / parallel copy) ────────────────────────────


class TestMergeArchiveUncompressed:

    def test_merge_uncompressed(self, workdir):
        kg = _make_kg(2)
        kg.archive("ds_unc", compress=False)

        merged = KnowledgeGraph()
        merged.merge_archive("ds_unc", compress=False)
        assert merged.n_samples == 2

    def test_parallel_copy_files_present(self, workdir):
        kg = _make_kg(3)
        kg.archive("ds_unc_par", compress=False)

        merged = KnowledgeGraph()
        merged.merge_archive("ds_unc_par", compress=False)

        ss_files = _structure_store_files(merged)
        path_basenames = _has_path_basenames(merged)
        for bn in path_basenames:
            assert bn in ss_files


# ── edge cases ──────────────────────────────────────────────────────────────


class TestArchiveEdgeCases:

    def test_archive_single_structure(self, workdir):
        kg = _make_kg(1)
        kg.archive("single")
        kg2 = KnowledgeGraph.unarchive("single.tar.gz")
        assert kg2.n_samples == 1

    def test_archive_already_exists_raises(self, workdir):
        kg = _make_kg(1)
        os.mkdir("existing_pkg")
        with pytest.raises(ValueError, match="already exists"):
            kg.archive("existing_pkg")

    def test_archive_tarball_already_exists_raises(self, workdir):
        kg = _make_kg(1)
        with open("existing_pkg.tar.gz", "w") as f:
            f.write("dummy")
        with pytest.raises(ValueError, match="tarball already exists"):
            kg.archive("existing_pkg")
