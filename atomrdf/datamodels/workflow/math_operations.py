"""
ASMO math-operation datamodels.
"""

from __future__ import annotations

from typing import List, Optional, Union
import uuid

from pydantic import Field
from rdflib import URIRef, Literal, XSD, RDF

from atomrdf.datamodels.activity import Activity
from atomrdf.datamodels.workflow.property import CalculatedProperty
from atomrdf.namespace import ASMO, PROV
from atomrdf.utils import get_sample_object


# ---------------------------------------------------------------------- #
# Internal helpers                                                         #
# ---------------------------------------------------------------------- #


def _resolve_operand(
    operand: Union[str, float],
    property_map: Optional[dict] = None,
) -> Union[URIRef, Literal]:
    """Convert a YAML operand value to an RDF term.

    Parameters
    ----------
    operand:
        Either a property-ID string (looked up in *property_map*) or a
        numeric scalar.
    property_map:
        Optional mapping of YAML local IDs to KG URI strings.
    """
    if isinstance(operand, str):
        if property_map:
            uri_str = property_map.get(operand, operand)
        else:
            uri_str = operand
        return URIRef(uri_str)
    return Literal(float(operand), datatype=XSD.float)


def _operand_from_node(node) -> Optional[Union[str, float]]:
    """Inverse of :func:`_resolve_operand`.

    Returns ``float`` for an ``xsd:float`` Literal; the URI string for a
    URIRef; ``None`` if *node* is ``None``.
    """
    if node is None:
        return None
    if isinstance(node, Literal):
        return float(node)
    return str(node)


def _add_result_to_graph(
    graph,
    activity: URIRef,
    result: CalculatedProperty,
) -> URIRef:
    result_uri = result.to_graph(graph)  # sets result.id, returns URIRef
    graph.add((result_uri, PROV.wasGeneratedBy, activity))
    graph.add((result_uri, ASMO.wasCalculatedBy, activity))
    if result.associate_to_sample:
        for sample_id in result.associate_to_sample:
            graph.add((URIRef(sample_id), ASMO.hasCalculatedProperty, result_uri))
    return result_uri


class Subtraction(Activity):
    """ASMO Subtraction: difference = minuend − subtrahend."""

    minuend: Optional[Union[str, float]] = Field(
        default=None,
        description="Minuend – property ID string or numeric scalar",
    )
    subtrahend: Optional[Union[str, float]] = Field(
        default=None,
        description="Subtrahend – property ID string or numeric scalar",
    )
    result: Optional[CalculatedProperty] = Field(
        default=None, description="Resulting calculated property"
    )

    def to_graph(self, graph, property_map: Optional[dict] = None) -> str:
        activity_id = f"subtraction:{uuid.uuid4()}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.Subtraction)

        if self.minuend is not None:
            graph.add(
                (
                    activity,
                    ASMO.hasMinuend,
                    _resolve_operand(self.minuend, property_map),
                )
            )
        if self.subtrahend is not None:
            graph.add(
                (
                    activity,
                    ASMO.hasSubtrahend,
                    _resolve_operand(self.subtrahend, property_map),
                )
            )
        if self.result is not None:
            _add_result_to_graph(graph, activity, self.result)

        return activity_id

    @classmethod
    def from_graph(cls, graph, activity_id: str) -> "Subtraction":
        activity = get_sample_object(activity_id)

        minuend_node = graph.value(activity, ASMO.hasMinuend)
        subtrahend_node = graph.value(activity, ASMO.hasSubtrahend)

        result_uri = next(
            (s for s, _, _ in graph.triples((None, PROV.wasGeneratedBy, activity))),
            None,
        )
        result = (
            CalculatedProperty.from_graph(graph, result_uri) if result_uri else None
        )

        return cls(
            id=str(activity_id),
            minuend=_operand_from_node(minuend_node),
            subtrahend=_operand_from_node(subtrahend_node),
            result=result,
        )


class Addition(Activity):
    """ASMO Addition: sum = addend₁ + addend₂ + …"""

    addend: Optional[List[Union[str, float]]] = Field(
        default=None,
        description="List of addends – each a property ID string or numeric scalar",
    )
    result: Optional[CalculatedProperty] = Field(
        default=None, description="Resulting calculated property"
    )

    def to_graph(self, graph, property_map: Optional[dict] = None) -> str:
        activity_id = f"addition:{uuid.uuid4()}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.Addition)

        for item in self.addend or []:
            graph.add((activity, ASMO.hasAddend, _resolve_operand(item, property_map)))
        if self.result is not None:
            _add_result_to_graph(graph, activity, self.result)

        return activity_id

    @classmethod
    def from_graph(cls, graph, activity_id: str) -> "Addition":
        activity = get_sample_object(activity_id)

        addends = [
            _operand_from_node(n) for n in graph.objects(activity, ASMO.hasAddend)
        ]

        result_uri = next(
            (s for s, _, _ in graph.triples((None, PROV.wasGeneratedBy, activity))),
            None,
        )
        result = (
            CalculatedProperty.from_graph(graph, result_uri) if result_uri else None
        )

        return cls(
            id=str(activity_id),
            addend=addends or None,
            result=result,
        )


class Multiplication(Activity):
    """ASMO Multiplication: product = factor₁ × factor₂ × …"""

    factor: Optional[List[Union[str, float]]] = Field(
        default=None,
        description="List of factors – each a property ID string or numeric scalar",
    )
    result: Optional[CalculatedProperty] = Field(
        default=None, description="Resulting calculated property"
    )

    def to_graph(self, graph, property_map: Optional[dict] = None) -> str:
        activity_id = f"multiplication:{uuid.uuid4()}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.Multiplication)

        for item in self.factor or []:
            graph.add((activity, ASMO.hasFactor, _resolve_operand(item, property_map)))
        if self.result is not None:
            _add_result_to_graph(graph, activity, self.result)

        return activity_id

    @classmethod
    def from_graph(cls, graph, activity_id: str) -> "Multiplication":
        activity = get_sample_object(activity_id)

        factors = [
            _operand_from_node(n) for n in graph.objects(activity, ASMO.hasFactor)
        ]

        result_uri = next(
            (s for s, _, _ in graph.triples((None, PROV.wasGeneratedBy, activity))),
            None,
        )
        result = (
            CalculatedProperty.from_graph(graph, result_uri) if result_uri else None
        )

        return cls(
            id=str(activity_id),
            factor=factors or None,
            result=result,
        )


class Division(Activity):
    """ASMO Division: quotient = dividend ÷ divisor."""

    dividend: Optional[Union[str, float]] = Field(
        default=None,
        description="Dividend – property ID string or numeric scalar",
    )
    divisor: Optional[Union[str, float]] = Field(
        default=None,
        description="Divisor – property ID string or numeric scalar",
    )
    result: Optional[CalculatedProperty] = Field(
        default=None, description="Resulting calculated property"
    )

    def to_graph(self, graph, property_map: Optional[dict] = None) -> str:
        activity_id = f"division:{uuid.uuid4()}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.Division)

        if self.dividend is not None:
            graph.add(
                (
                    activity,
                    ASMO.hasDividend,
                    _resolve_operand(self.dividend, property_map),
                )
            )
        if self.divisor is not None:
            graph.add(
                (
                    activity,
                    ASMO.hasDivisor,
                    _resolve_operand(self.divisor, property_map),
                )
            )
        if self.result is not None:
            _add_result_to_graph(graph, activity, self.result)

        return activity_id

    @classmethod
    def from_graph(cls, graph, activity_id: str) -> "Division":
        activity = get_sample_object(activity_id)

        dividend_node = graph.value(activity, ASMO.hasDividend)
        divisor_node = graph.value(activity, ASMO.hasDivisor)

        result_uri = next(
            (s for s, _, _ in graph.triples((None, PROV.wasGeneratedBy, activity))),
            None,
        )
        result = (
            CalculatedProperty.from_graph(graph, result_uri) if result_uri else None
        )

        return cls(
            id=str(activity_id),
            dividend=_operand_from_node(dividend_node),
            divisor=_operand_from_node(divisor_node),
            result=result,
        )


class Exponentiation(Activity):
    """ASMO Exponentiation: result = base ** exponent."""

    base: Optional[Union[str, float]] = Field(
        default=None,
        description="Base – property ID string or numeric scalar",
    )
    exponent: Optional[Union[str, float]] = Field(
        default=None,
        description="Exponent – property ID string or numeric scalar",
    )
    result: Optional[CalculatedProperty] = Field(
        default=None, description="Resulting calculated property"
    )

    def to_graph(self, graph, property_map: Optional[dict] = None) -> str:
        activity_id = f"exponentiation:{uuid.uuid4()}"
        self.id = activity_id
        activity = graph.create_node(activity_id, ASMO.Exponentiation)

        if self.base is not None:
            graph.add(
                (activity, ASMO.hasBase, _resolve_operand(self.base, property_map))
            )
        if self.exponent is not None:
            graph.add(
                (
                    activity,
                    ASMO.hasExponent,
                    _resolve_operand(self.exponent, property_map),
                )
            )
        if self.result is not None:
            _add_result_to_graph(graph, activity, self.result)

        return activity_id

    @classmethod
    def from_graph(cls, graph, activity_id: str) -> "Exponentiation":
        activity = get_sample_object(activity_id)

        base_node = graph.value(activity, ASMO.hasBase)
        exponent_node = graph.value(activity, ASMO.hasExponent)

        result_uri = next(
            (s for s, _, _ in graph.triples((None, PROV.wasGeneratedBy, activity))),
            None,
        )
        result = (
            CalculatedProperty.from_graph(graph, result_uri) if result_uri else None
        )

        return cls(
            id=str(activity_id),
            base=_operand_from_node(base_node),
            exponent=_operand_from_node(exponent_node),
            result=result,
        )


MATH_OPERATION_MAP = {
    "Subtraction": Subtraction,
    "Addition": Addition,
    "Multiplication": Multiplication,
    "Division": Division,
    "Exponentiation": Exponentiation,
}
