import enum
import bisect
import itertools
from collections import deque
from typing import Union, List


class SuffixTreeBuildingMethod(enum.Enum):
    naive = enum.auto()
    ukkonen = enum.auto()


class _TreeElem:
    def __init__(self, data=None, index=None, children_type=list):
        self.node = data
        if children_type == list:
            self.children = [_TreeElem(index)] if index else []
        elif children_type == dict:
            self.children = {index: _TreeElem(index)} if index else {}

    def add(self, suffix):
        current = self
        while True:
            pos = bisect.bisect(current.children, suffix)
            if not self._is_prefix_presented(current, pos, suffix):
                pos = pos - 1
            if self._is_prefix_presented(current, pos, suffix):
                general = current.children[pos].extract_general_prefix(suffix)
                if general < len(current.children[pos].node):  # added suffix and current vertex has common prefix
                    self._split_and_insert_common_prefix(current, general, pos)
                if general < len(suffix.node):  # vertex is a prefix case
                    suffix.node = suffix.node[general:]
                else:
                    # vertex is already presented, we have to add just index
                    bisect.insort(current.children[pos].children, suffix.children[0])
                    break
                current = current.children[pos]
            else:
                bisect.insort(current.children, suffix)
                break

    def _is_prefix_presented(self, current, pos, suffix):
        return pos > -1 and pos < len(current.children) \
            and type(current.children[pos].node) == str and suffix.node.startswith(current.children[pos].node[0])

    def _split_and_insert_common_prefix(self, current, general, pos):
        rest = _TreeElem(current.children[pos].node[general:])
        rest.children = current.children[pos].children
        current.children[pos].node = current.children[pos].node[:general]
        current.children[pos].children = [rest]

    def __str__(self):
        children = ", ".join(str(child) if type(child) == list else str(self.children[child])
                             for child in self.children)
        return f"Elem(node={self.node}, children={children})"

    def __eq__(self, other):
        return self.node == other.node

    def __lt__(self, other):
        return (type(self.node) == int, self.node) < (type(other.node) == int, other.node)

    def __gt__(self, other):
        return (type(self.node) == int, self.node) > (type(other.node) == int, other.node)

    def extract_general_prefix(self, other):
        for general, (a, b) in enumerate(zip(self.node, other.node)):
            if a != b:
                break
        else:
            general += 1
        return general


class SuffixTree:
    def __init__(self, collection: Union[str, List[str]], method=SuffixTreeBuildingMethod.ukkonen):
        self._n_string = 0
        self.method = method
        if method == SuffixTreeBuildingMethod.naive:
            self._naive(collection)
        elif method == SuffixTreeBuildingMethod.ukkonen:
            self._ukkonen(collection)
        else:
            raise ValueError("Wrong method, please choose one from SuffixTreeBuildingMethod")

    def _naive(self, collection):
        collection = [collection] if isinstance(collection, str) else collection
        self._tree = _TreeElem()

        for elem in collection:
            self.add(elem)

    def _ukkonen(self, collection):
        self._n_string += 1
        collection = [collection] if isinstance(collection, str) else collection
        self._tree = _TreeElem(children_type=dict)
        self._first_elem = collection[0]

        for elem in collection:
            self._add_ukkonen(elem)

    def _add_ukkonen(self, elem):
        active_edge = self._tree
        active_length = 0
        for index, symbol in enumerate(itertools.chain(elem, (self._n_string,))):
            if not active_length:
                if symbol not in active_edge.children:
                    child = _TreeElem((index, None), children_type=dict, index=self._n_string)
                    active_edge.children[symbol] = child
                else:
                    active_node = active_edge.children[symbol].node[0]
                    active_edge = active_edge.children[symbol]
                    active_length += 1
            else:
                # todo: check if current suffix is longer then current node
                print(active_node, active_length)
                if active_node + active_length < len(self._first_elem)\
                        and self._first_elem[active_node + active_length] != symbol:
                    while active_length:
                        self._split_node_ukkonen(active_edge, active_length, active_node, elem, index, symbol)

                        active_length -= 1
                        active_node += 1
                        active_edge = self._tree.children[elem[active_node]]
                else:
                    active_length += 1

    def _split_node_ukkonen(self, active_edge, active_length, active_node, elem, index, symbol):
        start, end = active_edge.node
        active_edge.node = start, start + active_length

        #assert not end or end >= active_node + active_length
        child = self._first_elem[active_node + active_length]
        child_tree = _TreeElem((start + active_length, end), children_type=dict)
        child_tree.children = active_edge.children
        active_edge.children = {child: child_tree}

        child = (index, None) if type(symbol) != int else symbol
        active_edge.children[symbol] = _TreeElem(child, children_type=dict)

    def add(self, string):
        self._n_string += 1
        cur = ""

        if self.method == SuffixTreeBuildingMethod.naive:
            for symbol in reversed(string):
                cur = symbol + cur
                self._tree.add(_TreeElem(cur, self._n_string))
        else:
            self._add_ukkonen(string)


    def make_common_tree(self):
        self._n_string = 1
        self._tree.children = [child for child in self._tree.children
                               if not (len(child.children) == 1 and type(child.children[0].node) == int)]
        stack = [(self._tree, child, child.node, False) for child in self._tree.children]
        longests = []
        while stack:
            parent, current, string, visited = stack[-1]
            stack[-1] = parent, current, string, True
            if len(current.children) == 1 and type(current.children[0].node) == int:
                self._remove_uncommon_vertex(current, parent)
            else:
                longests = self._current_longests(longests, string, visited)
            if visited or len(current.children) == 1:
                stack.pop()
                self._remove_endings(current)
            else:
                stack.extend((current, child, string + child.node, False)
                             for child in current.children if type(child.node) != int)
        return longests

    def _remove_uncommon_vertex(self, current, parent):
        parent.children.pop(parent.children.index(current))
        child = current.children[0]
        pos = bisect.bisect_left(parent.children, child)
        if not (0 <= pos < len(parent.children) and parent.children[pos] == child):
            bisect.insort(parent.children, child)

    def _current_longests(self, longests, string, visited):
        if visited:
            longest_len = len(longests[0]) if longests else 0
            if len(string) > longest_len:
                longests = [string]
            elif len(string) == longest_len:
                longests.append(string)
        return longests

    def _remove_endings(self, current):
        for index in range(-1, -len(current.children), -1):
            node = current.children[index].node
            if type(node) == int and node != 1:
                current.children.pop()

    def to_graphviz(self):
        number = 1
        deq = deque([(self._tree, number)])
        nodes = ""
        edges = ""
        while deq:
            current, parent_number = deq.popleft()
            node_label = current.node if type(current.node) != tuple else (current.node, self._first_elem[current.node[0]:current.node[1]])
            nodes += str(parent_number) + f' [label="{node_label}"]\n'
            iterable = current.children.values() if type(current.children) == dict else current.children
            for current_index, child in enumerate(iterable, 1):
                number += 1
                deq.append((child, number))
                edges += str(parent_number) + "--" + str(number) + "\n"
        return "strict graph G {\n" + nodes + edges + "}"


    def __str__(self):
        return str(self._tree)

    def __repr__(self):
        return repr(self._tree)
