import matplotlib.pyplot as plt
import numpy as np
import argparse
import statistics as stat
import re
import sys
import os
import logging

from ete3 import Tree
from typing import List
from collections import namedtuple
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.ensemble import RandomForestClassifier

# Removing Outliers Per Rank
RED_RANK = namedtuple('RedRank', 'red rank')


def clean_lineage_string(lineage: str):
    """
    Removes superfluous taxonomic ranks and characters that make lineage comparisons difficult

    :param lineage: A taxonomic lineage string where each rank is separated by a semi-colon
    :return: String with the purified taxonomic lineage adhering to the NCBI hierarchy
    """
    non_standard_names_re = re.compile(" group| cluster| complex", re.IGNORECASE)
    bad_strings = ["cellular organisms; ", "delta/epsilon subdivisions; ", "\(miscellaneous\)", "[a-p]__"]
    for bs in bad_strings:
        lineage = re.sub(bs, '', lineage)
    # filter 'group' and 'cluster'
    if non_standard_names_re.search(lineage):
        reconstructed_lineage = ""
        ranks = lineage.split("; ")
        for rank in ranks:
            if not non_standard_names_re.search(rank):
                reconstructed_lineage = reconstructed_lineage + str(rank) + '; '
        reconstructed_lineage = re.sub('; $', '', reconstructed_lineage)
        lineage = reconstructed_lineage
    return lineage


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
        cleaned_lineage = clean_lineage_string(lineage.rstrip())
        tax_path = cleaned_lineage.split('; ')
        tree_num_to_lineage_map[tree_id] = tax_path
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
        node.add_features(_pass=True)

    return t


# Functions to find avg dist from leaves to a nodes parent (for calculation of RED)
class Dist(object):
    @staticmethod
    def avg_dist_to_this_node(node):
        """
        get avg distance from leaves to this node's parent
        """
        return Dist.avg(Dist.get_list_distances_of_this_nodes_leaves_to_node(node))

    @staticmethod
    def avg(l: List[float]):
        """
        averages list of numbers
        """
        if len(l) == 0:
            return None
        else:
            return sum(l)/len(l)

    @staticmethod
    def get_list_distances_of_this_nodes_leaves_to_node(node):
        """
        from this node, get list of distances from child leaves
        """
        distances = []
        for n in node.get_leaves():
            if n._pass:
                distances.append(n.get_distance(node))
        return distances


# Functions for getting LCA from leaf node info
class LCA(object):
    skip_taxa = ["cellular organisms"]
    @staticmethod
    def get_lca_lineage(node):
        """
        main function that returns lineage based on 'passed' leaf nodes
        """
        return LCA.get_commons_from_mult_lists(LCA.get_leaf_linages_passed(node))

    @staticmethod
    def assign_pass(t, keyword_parameters):
        for leaf in t:
            if 'remove_strings' in keyword_parameters:
                for word in keyword_parameters['remove_strings']:
                    if leaf in LCA.list_leaves_at_rank(t, word, 10):
                        leaf.add_features(_pass=False)
                        logging.debug("removed = " + str(leaf.lineage) + "\n")
            if "min_lin_depth" in keyword_parameters and keyword_parameters["min_lin_depth"] > 0:
                min_lin_depth = keyword_parameters["min_lin_depth"]
                if LCA.lineage_length(leaf.lineage) < min_lin_depth:
                        leaf.add_features(_pass=False)
                        logging.debug("removed = " + str(leaf.lineage) + "\n")
        return 0

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
            if LCA.lineage_length(node.lineage) <= rank:
                for tax in node.lineage:
                    if cls in tax:
                        l.append(node)
        return l

    @staticmethod
    def lineage_length(lineage_list) -> int:
        """
        assesses rank [1:8] where 1 is the domain/kingdom, 2 is the phylum... 7 being species
        """
        if lineage_list is None:
            return 0
        else:
            d = 0
            for taxon in lineage_list:
                if taxon not in LCA.skip_taxa:
                    d += 1
            return d

    @staticmethod
    def get_commons_from_mult_lists(l: list):
        """
        returns a list of common elements from multiple strings
        """
        if l is None:
            return []
        elif type(l) is list and len(l) > 0:
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
            elif node.is_leaf():
                node.add_features(red=1)
        for leaf in t:
            if LCA.lineage_length(leaf.lineage) >= 7:
                leaf.add_features(red=0.9999)
        return t

    @staticmethod
    def label_red(node):
        """labels RED to a single node given parent has RED value"""
        return node.add_features(red=RED.get_red(node))

    @staticmethod
    def get_red(node):
        """gets the RED value associated with this node given parent has RED value"""
        if Dist.avg_dist_to_this_node(node) is None:
            return None
        else:
            a = node.get_distance(node.up)
            b = Dist.avg_dist_to_this_node(node)
            x = node.up.red
            if a + b != 0:
                red = x + (a/(a + b))*(1 - x)
            else:
                red = x
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
        if not l:
            logging.debug("No descendants of node '" + str(t.name) + "'.\n")
            return 0
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
                if node.lineage:
                    node.add_features(rank=LCA.lineage_length(node.lineage))
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
        logging.info("\n".join([str(k) + " = '" + str(v) + "'" for k, v in kwargs_dict.items()]) + "\n")
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


def cull_outliers(data: list, dev=3):
    """
    Returns the Interquartile Range (IQR) of a list after filtering outliers
    based on log transformed data where outliers are farther than 1 std-dev from the median and
    an un-transformed distribution where outliers are farther than `dev` standard deviations from the median.

    :param data: A list of floats
    :param dev: Number of acceptable deviations from the median; beyond this, values are outliers and removed
    :return: A smaller list of floats
    """
    # Reject outliers from ln-transformed distribution
    ln_a = np.log10(1.0 * np.array(data))
    noo_a = np.power(10, ln_a[abs(ln_a - np.median(ln_a)) < 2 * np.std(ln_a)])

    # Reject outliers from untransformed distribution
    d = np.abs(noo_a - np.median(noo_a))
    mdev = np.median(d)
    s = d / mdev if mdev else 0
    noo_a = noo_a[s < dev]
    try:
        if isinstance(noo_a[0], np.ndarray):
            return list(noo_a[0])
        else:
            return list(noo_a)
    except IndexError:
        return []


def list_nodes_of_rank(t, rank):
    """
    Generates a list of nodes of a given rank integer.

    :param t: tree
    :param rank: integer, range[1:7]
    :return: list
    """
    nodes_of_rank = []
    for node in t.iter_descendants():
        if node.rank == rank:
            nodes_of_rank.append(node)
    return nodes_of_rank


def dict_of_nodes_of_rank(t: Tree, rank):
    """
    Generates dictionary of nodes included in a rank integer. RED values and Ranks
    become the value of their respective key (node) in namedtuple format: RED_RANK(RED, rank).
    :param t: tree
    :param rank: int
    :return: dictionary
    """
    nodes_of_rank = list_nodes_of_rank(t, rank)
    dict = {}
    for node in nodes_of_rank:
        if node.red is not None:
            if node.red < 1:
                if node.rank is not None:
                    RR = RED_RANK(node.red, node.rank)
                    dict[node] = RR
    return dict


def list_inliers_outliers(nodes: dict):
    """
    Returns dictionaries of inlier nodes, lower outlier nodes, upper outlier nodes from
    a dictionary of nodes (here we use the inputs as dictionaries of nodes from a single rank)

    :param nodes: dictionary
    :return: three dictionaries
    """
    reds = []
    for RR in nodes.values():
        reds.append(RR.red)
    noo_a = cull_outliers(reds)
    median = np.median(noo_a)

    inliers = {}
    low_outliers = {}
    high_outliers = {}
    for node in nodes:
        if nodes[node].red in noo_a:
            inliers[node] = RED_RANK(node.red, node.rank)
        else:
            if nodes[node].red < median:
                low_outliers[node] = RED_RANK(node.red, node.rank)
            elif nodes[node].red > median:
                high_outliers[node] = RED_RANK(node.red, node.rank)

    return inliers, low_outliers, high_outliers


def get_full_inliers_and_outliers(t, r1, r2):
    """
    Returns two dictionaries of nodes that are either inliers or outliers relative to the
    red values of nodes in the same rank.

    :param t: tree
    :param r1: int, bottom rank range to include
    :param r2: int, upper rank range to include
    :return: two dictionaries of nodes as keys, RED_RANK(red, rank) as values
    """
    full_outliers = {}
    full_inliers = {}
    for num in range(r1, r2+1):
        rank_dict = dict_of_nodes_of_rank(t, num)
        inliers, low_outliers, high_outliers = list_inliers_outliers(rank_dict)
        for node in low_outliers:
            full_outliers[node] = RED_RANK(node.red, node.rank)
        for node in high_outliers:
            full_outliers[node] = RED_RANK(node.red, node.rank)
        for node in inliers:
            full_inliers[node] = RED_RANK(node.red, node.rank)
    return full_inliers, full_outliers


def plot_fit(model, x_train, x_test, y_train, y_test):
    fit = model.fit(x_train, y_train)
    model.score(x_test, y_test)
    lex = []
    ley = []
    for num in range(0, 100):
        lex.append(num / 100)
        ley.append(fit.predict(np.array(num / 100).reshape(-1, 1))[0])
    plt.plot(lex, ley, 'b-')
    return


def model_graph(t, model_type):
    """
    Graphically fits a model using the ML to processed dataset of a tree. Model types include
    linear (linear regression), logistic (logistic regression), and forest (random forest
    classifier).

    :param t: tree
    :param model_type: str
    :return: matplotlib graph
    """
    inliers, outliers = get_full_inliers_and_outliers(t, 1, 7)
    reds = []
    ranks = []
    for node in inliers:
        reds.append(inliers[node].red)
        ranks.append(inliers[node].rank)
    x_train, x_test, y_train, y_test = train_test_split(np.array(reds).reshape(-1, 1), ranks,
                                                        test_size=0.2, random_state=0)

    if model_type == 'random_forest':
        model = RandomForestClassifier(class_weight='balanced', min_weight_fraction_leaf=0.25, n_estimators=100)
    elif model_type == 'linear':
        model = LinearRegression()
    elif model_type == 'logistic':
        model = LogisticRegression(multi_class='multinomial', solver='lbfgs', class_weight='balanced')
    else:
        logging.error("Unknown model type '" + str(model_type) + "'.\n")
        sys.exit(3)

    plot_fit(model, x_train, x_test, y_train, y_test)
    score = model.score(x_test, y_test)
    model_score = str(score)
    logging.info("Model score = " + model_score + "\n")
    plt.text(x=0.1, y=6.5, s=0, text='model score = ' + model_score)
    return plt.show()


# Graphing red vs rank
def graph_red_vs_rank(t, model_type):
    """returns red vs rank graphic with median values at each rank and the correlation coefficient"""
    inliers, outliers = get_full_inliers_and_outliers(t, 1, 7)
    reds = []
    ranks = []
    for node in inliers:
        reds.append(inliers[node].red)
        ranks.append(inliers[node].rank)
    plt.xlabel('RED')
    plt.ylabel('Rank')
    plt.title('RED Assignment')
    plt.plot(reds, ranks, 'ro')
    plt.axis([0, 1, 0, 8])
    for num in range(0, 8):
        if RED.median_red(t, num) is None:
            continue
        else:
            median = RED.median_red(t, num)
            plt.plot(median, num, 'g^')
            plt.text(x=median-0.05, y=num+0.4, s=0, text=str(round(median, 4)))
    plt.grid(True)
    model_graph(t, model_type)
    return plt.show()


def get_arguments():
    parser = argparse.ArgumentParser(description='Calculate average distance of taxonomic rank of a tree to the root')
    parser.add_argument('-i', '--tax_ids',
                        type=str, metavar='tax_ids.txt', required=True,
                        help="Taxonomic IDs of a tree file")
    parser.add_argument('-t', '--tree',
                        type=str, metavar='tree.txt', required=True,
                        help="Tree file of interest, in Newick format")
    parser.add_argument('-m', '--model',
                        type=str, metavar='', required=False, default="linear",
                        choices=["linear", "logistic", "random_forest"],
                        help="Model to fit to the RED distances across taxonomic ranks [ DEFAULT = linear ]")
    parser.add_argument('-r', '--remove',
                        type=str, metavar='', required=False, nargs='+', default="",
                        help="Remove lineages containing specific strings "
                             "(e.g. metagenome, unclassified, candidatus, environmental) at given rank")
    parser.add_argument('-l', '--lineage_len',
                        type=int, metavar='', required=False, default=2,
                        help="The minimum taxonomic lineage resolution/depth (Kingdom = 1, Phylum = 2, etc.). "
                             "Removes leaves with a truncated lineage at this depth [ DEFAULT = 2 ]")
    args = parser.parse_args()
    return args


class MyFormatter(logging.Formatter):

    error_fmt = "%(levelname)s - %(module)s, line %(lineno)d:\n%(message)s"
    warning_fmt = "%(levelname)s:\n%(message)s"
    debug_fmt = "%(asctime)s\n%(message)s"
    info_fmt = "%(message)s"

    def __init__(self):
        super().__init__(fmt="%(levelname)s: %(message)s",
                         datefmt="%d/%m %H:%M:%S")

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._style._fmt = MyFormatter.debug_fmt

        elif record.levelno == logging.INFO:
            self._style._fmt = MyFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._style._fmt = MyFormatter.error_fmt

        elif record.levelno == logging.WARNING:
            self._style._fmt = MyFormatter.warning_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result


def prep_logging(log_file_name=None, verbosity=False):
    """
    Allows for multiple file handlers to be added to the root logger, but only a single stream handler.
    The new file handlers must be removed outside of this function explicitly
    :param log_file_name:
    :param verbosity:
    :return:
    """
    if verbosity:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO

    # Detect whether a handlers are already present and return if true
    logger = logging.getLogger()
    if len(logger.handlers):
        return

    formatter = MyFormatter()
    # Set the console handler normally writing to stdout/stderr
    ch = logging.StreamHandler()
    ch.setLevel(logging_level)
    ch.terminator = ''
    ch.setFormatter(formatter)

    if log_file_name:
        output_dir = os.path.dirname(log_file_name)
        try:
            if output_dir and not os.path.isdir(output_dir):
                os.makedirs(output_dir)
        except (IOError, OSError):
            sys.stderr.write("ERROR: Unable to make directory '" + output_dir + "'.\n")
            sys.exit(3)
        logging.basicConfig(level=logging.DEBUG,
                            filename=log_file_name,
                            filemode='w',
                            datefmt="%d/%m %H:%M:%S",
                            format="%(asctime)s %(levelname)s:\n%(message)s")
        logging.getLogger('').addHandler(ch)
        logging.getLogger('').propagate = False
    else:
        logging.basicConfig(level=logging_level,
                            datefmt="%d/%m %H:%M:%S",
                            format="%(asctime)s %(levelname)s:\n%(message)s")
    return


if __name__ == '__main__':
    args = get_arguments()
    prep_logging("./RedAssignment_log.txt")
    ref_tree = open_tree_file(args.tax_ids, args.tree)
    Map.class_all_nodes(ref_tree,
                        min_lin_depth=args.lineage_len,
                        remove_strings=args.remove)
    RED.apply_all(ref_tree)
    Map.label_nodes(ref_tree)
    graph_red_vs_rank(ref_tree, args.model)
