import pyham
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")

nwk_path = "Output/EstimatedSpeciesTree.nwk"
tree_str = pyham.utils.get_newick_string(nwk_path, type="nwk")
orthoxml_path = "Output/HierarchicalGroups.orthoxml"

ham_analysis = pyham.Ham(nwk_path, orthoxml_path, use_internal_name=False)

treeprofile = ham_analysis.create_tree_profile(outfile="tp.html")

