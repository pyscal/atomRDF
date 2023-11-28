def create_graph():
    n = 5
    from pyscal_rdf import StructureGraph
    g = StructureGraph(store='SQLAlchemy', store_file="testfile_multi_2.db")
    #g = StructureGraph()
    struct_Fe = g.create_structure("l12", element=['Al', 'Ni'], lattice_constant=3.57, repetitions=[n,n,n])
    g.add_structure_to_graph(struct_Fe)
    g.add_structure_to_graph(struct_Fe)
    g.add_structure_to_graph(struct_Fe)
