import matplotlib.pyplot as plt
import numpy as np
import argparse
import statistics as stat

from ete3 import Tree
from typing import List


# Opening Tree and ID Files
def open_tree_file(taxa_table_file: str, newick_file: str):
    """
    There is some sort of taxonomic lineage-based filtering
    :param taxa_table_file:
    :param newick_file:
    :return:
    """
    taxa_table_file = str(taxa_table_file)
    newick_file = str(newick_file)

    tax_ids_handler = open(taxa_table_file)
    tree_num_to_lineage_map = {}
    for line in tax_ids_handler:
        tree_id, desc, lineage = line.split('\t')
        tax = lineage.rstrip().split('; ')
        tree_num_to_lineage_map[tree_id] = tax
    tax_ids_handler.close()

    tree_handler = open(newick_file)
    tree_contents = tree_handler.read()
    t = Tree(tree_contents)
    tree_handler.close()

    # Mapping lineages to leaf nodes
    for node in t.get_leaves():
        num = node.name
        node.add_features(lineage=tree_num_to_lineage_map[str(num)])

    # Add cellular organisms tax group to lineages w/o it
    for node in t.traverse():
        if node.is_leaf() and 'cellular organisms' not in node.lineage:
            new_lin = ['cellular organisms'] + node.lineage
            node.add_features(lineage=new_lin)
            # print(node.lineage)
    return t


# Functions to find avg dist from leaves to a nodes parent (for calculation of RED)
class Dist(object):
    @staticmethod
    def avg_dist_to_this_nodes_parent(node):
        """
        get avg distance from leaves to this node's parent
        """
        return Dist.avg(Dist.get_list_distances_of_this_nodes_leaves_to_nodes_parent(node))

    @staticmethod
    def avg(l: List[float]):
        """
        averages list of numbers
        """
        return sum(l)/len(l)

    @staticmethod
    def get_list_distances_of_this_nodes_leaves_to_nodes_parent(node):
        """
        from this node, get list of distances from child leaves
        """
        l = []
        parent_node = node.up
        for n in node.get_leaves():
            l.append(n.get_distance(parent_node))
        return l


# Functions for getting LCA from leaf node info
class LCA(object):

    @staticmethod
    def get_lca_lineage(node):
        """
        main function that returns lineage based on 'passed' leaf nodes
        """
        return LCA.get_commons_from_mult_lists(LCA.get_leaf_linages_passed(node))

    @staticmethod
    def assign_pass(t, keyword_parameters):
        for node in t.traverse():
            node.add_features(_pass=True)
        for leaf in t:
            if "min_lin_depth" in keyword_parameters and keyword_parameters["min_lin_depth"] > 0:
                if LCA.rank_assessment(leaf.lineage) < keyword_parameters["min_lin_depth"]:
                    # print("Too short:", leaf.lineage)
                    leaf.add_features(_pass=False)
            # TODO: Use the new 'remove_strings' argument, controlled by args.remove to remove all unwanted lineages
            if "remove_strings" in keyword_parameters and len(keyword_parameters["remove_strings"]) > 0:
                print(keyword_parameters["remove_strings"])
            # TODO: Remove all of these obsolete filters
            if 'r1leaves' in keyword_parameters:
                if keyword_parameters['r1leaves']:
                    if LCA.rank_assessment(leaf.lineage) == 1:
                        leaf.add_features(_pass=False)
            if 'r2leaves' in keyword_parameters:
                if LCA.rank_assessment(leaf.lineage) == 2:
                    if keyword_parameters['r2leaves']:
                        leaf.add_features(_pass=False)
            if 'metagenome' in keyword_parameters:
                if leaf in LCA.list_leaves_at_rank(t, 'metagenome', keyword_parameters['metagenome']):
                    leaf.add_features(_pass=False)
            if 'environmental_samples' in keyword_parameters:
                if leaf in LCA.list_leaves_at_rank(t, 'environmental samples', keyword_parameters['environmental_samples']):
                    leaf.add_features(_pass=False)
            if 'unclassified' in keyword_parameters:
                if leaf in LCA.list_leaves_at_rank(t, 'unclassified', keyword_parameters['unclassified']):
                    leaf.add_features(_pass=False)
            if 'candidatus' in keyword_parameters:
                if leaf in LCA.list_leaves_at_rank(t, 'Candidatus', keyword_parameters['candidatus']):
                    leaf.add_features(_pass=False)
            if 'miscellaneous' in keyword_parameters:
                if leaf in LCA.list_leaves_at_rank(t, 'miscellaneous', keyword_parameters['miscellaneous']):
                    leaf.add_features(_pass=False)
        return

    @staticmethod
    def get_leaf_linages_passed(node):
        """
        appends list of leaves that are 'passed'
        """
        l = []
        for n in node.get_leaves():
            if n._pass:
                l.append(n.lineage)
        return l

    @staticmethod
    def list_leaves_at_rank(t, cls, rank):
        """
        given tree and rank level you want to assess, returns a list of node names of a certain class
        (i.e. 'environmental samples') that is within that rank level and above
        e.g. t, cls='metagenome', rank=5 --> returns list of leaves at rank 1-5 that have metagenome within its lineage
        """
        l = []
        for node in t.get_leaves():
            if LCA.rank_assessment(node.lineage) <= rank:
                for tax in node.lineage:
                    if cls in tax:
                        l.append(node)
        return l

    @staticmethod
    def rank_assessment(l):
        """
        assesses rank [1:8] where 1 being cellular organism, 2 being domain... 8 being species
        """
        ni = 0
        if l is None or l == [] or l == 'None':  # Is the last condition necessary?
            return None
        else:
            for tax in l:
                if 'group' in tax or 'cluster' in tax:
                    ni = ni + 1
            rank = len(l) - ni
            if rank > 8:
                return 8
            else:
                return rank

    @staticmethod
    def get_commons_from_mult_lists(l):
        """
        returns a list of common elements from multiple strings
        """
        if l != [] or l is None:
            c = l[0]
            for num in range(1, (len(l))):
                c = LCA.get_commons_from_2_list(c, l[num])
            return c

    @staticmethod
    def get_commons_from_2_list(l1, l2):
        """
        returns a list of common elements from two lists
        """
        result = []
        for element in l1:
            if element in l2:
                result.append(element)
        return result


# RED functions and outputs
class RED(object):

    @staticmethod
    def apply_all(t):
        """applies RED to each node"""
        t.add_features(red=0)
        for node in t.iter_descendants('preorder'):
            if not node.is_leaf():
                RED.label_red(node)
                # print(node, node.red)
            elif node.is_leaf():
                node.add_features(red=1)
        return t

    @staticmethod
    def label_red(node):
        """labels RED to a single node given parent has RED value"""
        return node.add_features(red=RED.get_red(node))

    @staticmethod
    def get_red(node):
        """gets the RED value associated with this node given parent has RED value"""
        red = node.up.red + (node.get_distance(node.up)/Dist.avg_dist_to_this_nodes_parent(node))*(1 - node.up.red)
        return red

    @staticmethod
    def avg_red(t, rank):
        """
        returns average red value for a rank
        """
        l = []
        for node in t.iter_descendants():
            if node.rank == rank:
                if node.red < 1:
                    l.append(node.red)
        return Dist.avg(l)

    @staticmethod
    def median_red(t, rank):
        """
        returns median red value for rank
        """
        l = []
        for node in t.iter_descendants():
            if node.rank == rank:
                if node.red < 1:
                    if node.rank is not None:
                        l.append(node.red)
        if l == []:
            print(0)
        else:
            return stat.median(l)


# Adding lineages and ranks to internal nodes
class Map(object):

    @staticmethod
    def label_nodes(t):
        """
        label all node lineages based on rank level of filter
        and gives rank feature to all nodes
        note: nodes at level of filter only useful for distance measures
        """
        for node in t.traverse():
            try:
                if not node.lineage is None:
                    node.add_features(rank=LCA.rank_assessment(node.lineage))
                else:
                    node.add_features(rank=None)
            except AttributeError:
                node.add_features(rank=None)

    @staticmethod
    def class_all_nodes(t, **kwargs):
        """
        adds lineage feature to all nodes with available leaf descendants.
        LCA from leaf node info
        """
        kwargs_dict = kwargs
        LCA.assign_pass(t, kwargs_dict)
        Map.label_nodes(t)
        for node in reversed(list(t.traverse('levelorder'))):
            if node.is_leaf() is False:
                new_lin = LCA.get_lca_lineage(node)
                if new_lin is None:
                    node.add_features(lineage=None)
                else:
                    node.add_features(lineage=new_lin)
                    Map.label_nodes(t)


# Graphing red vs rank
def graph_red_vs_rank(t):
    """returns red vs rank graphic with median values at each rank and the correlation coefficient"""
    reds = []
    ranks = []
    for node in t.iter_descendants():
        if node.red < 1:
            if node.rank is not None:
                reds.append(node.red)
                ranks.append(node.rank)
    plt.xlabel('RED')
    plt.ylabel('Rank')
    plt.title('RED Assignment')
    plt.plot(reds, ranks, 'ro')
    plt.axis([0, 1, 1, 9])
    for num in range(2, 9):
        if RED.median_red(t, num) is None:
            continue
        else:
            median = RED.median_red(t, num)
            plt.plot(median, num, 'g^')
            plt.text(x=median-0.05, y=num+0.4, s=0, text=str(round(median, 4)))
    plt.grid(True)
    corrcoef = str(round(np.corrcoef(reds, ranks)[1, 0], 4))
    plt.text(x=0.1, y=6.5, s=0, text='corrcoef = ' + corrcoef)
    return plt.show()


def get_arguments():
    parser = argparse.ArgumentParser(description='Calculate average distance of taxonomic rank of a tree to the root')
    parser.add_argument('-i', '--tax_ids',
                        type=str, metavar='', required=True,
                        help="taxonomic IDs of a tree file")
    parser.add_argument('-t', '--tree',
                        type=str, metavar='', required=True,
                        help="tree file of interest")
    parser.add_argument('-r', '--remove',
                        type=str, metavar='', required=False, nargs='+', default="",
                        help="Remove lineages with strings "
                             "(e.g. metagenome, unclassified, candidatus, environmental) at given rank")
    parser.add_argument('-l', '--lineage_len',
                        type=int, metavar='', required=False, default=2,
                        help="The minimum taxonomic lineage resolution/depth (Kingdom = 1, Phylum = 2, etc.). "
                             "Removes leaves with a truncated lineage at this depth [ DEFAULT = 2 ]")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_arguments()
    t = open_tree_file(args.tax_ids, args.tree)
    RED.apply_all(t)
    Map.class_all_nodes(t,
                        min_lin_depth=args.lineage_len,
                        remove_strings=args.remove,
                        miscellaneous=8)
    Map.label_nodes(t)
    graph_red_vs_rank(t)
