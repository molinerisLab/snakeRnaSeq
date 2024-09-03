#!/usr/bin/env python
#
# Copyright 2008 Gabriele Sales <gbrsales@gmail.com>; 2014 Ivan Molineris <ivan.molineris@gmail.com>

from array import array
from optparse import OptionParser
from sys import stdin
from collections.abc import MutableSet
from collections import OrderedDict
from weakref import proxy

NA=None
ROW_ID=None

class Link(object):
    __slots__ = 'prev', 'next', 'key', '__weakref__'

#https://github.com/ActiveState/code/blob/3b27230f418b714bc9a0f897cb8ea189c3515e99/recipes/Python/576696_OrderedSet_with_Weakrefs/recipe-576696.py
class OrderedSet(MutableSet):
    'Set the remembers the order elements were added'
    # Big-O running times for all methods are the same as for regular sets.
    # The internal self.__map dictionary maps keys to links in a doubly linked list.
    # The circular doubly linked list starts and ends with a sentinel element.
    # The sentinel element never gets deleted (this simplifies the algorithm).
    # The prev/next links are weakref proxies (to prevent circular references).
    # Individual links are kept alive by the hard reference in self.__map.
    # Those hard references disappear when a key is deleted from an OrderedSet.

    def __init__(self, iterable=None):
        self.__root = root = Link()         # sentinel node for doubly linked list
        root.prev = root.next = root
        self.__map = {}                     # key --> link
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.__map)

    def __contains__(self, key):
        return key in self.__map

    def add(self, key):
        # Store new key in a new link at the end of the linked list
        if key not in self.__map:
            self.__map[key] = link = Link()            
            root = self.__root
            last = root.prev
            link.prev, link.next, link.key = last, root, key
            last.next = root.prev = proxy(link)

    def discard(self, key):
        # Remove an existing item using self.__map to find the link which is
        # then removed by updating the links in the predecessor and successors.        
        if key in self.__map:        
            link = self.__map.pop(key)
            link.prev.next = link.next
            link.next.prev = link.prev

    def __iter__(self):
        # Traverse the linked list in order.
        root = self.__root
        curr = root.next
        while curr is not root:
            yield curr.key
            curr = curr.next

    def __reversed__(self):
        # Traverse the linked list in reverse order.
        root = self.__root
        curr = root.prev
        while curr is not root:
            yield curr.key
            curr = curr.prev

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = next(reversed(self)) if last else next(iter(self))
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return not self.isdisjoint(other)


class SparseMatrix(object):
    def __init__(self, cols_item=OrderedSet()):
        self.data = OrderedDict() 
        self.cols_item = cols_item
        self.cols_item_sorted = []

    def set(self, x, y, v):
        pos = self.data.get(x,{})
        if y in pos:
            exit("only one value per position allowed, duplicated value for (%s,%s)" % (x,y))
        pos[y]=v
        self.data[x]=pos
        self.cols_item.add(y)

    def rows(self):
        for x in self.data.keys():
            retval = [x]
            for y in self.cols_item_sorted:
                retval.append(self.data[x].get(y,None))
            yield retval

    def head(self, sort=True):
        global ROW_ID
        if sort:
            self.cols_item_sorted = sorted(list(self.cols_item))
        else:
            self.cols_item_sorted = list(self.cols_item)
        return "%s\t%s" % (ROW_ID, "\t".join(self.cols_item_sorted))

    def clean(self):
        self.data={}

def None2NA(v):
    if v is None:
        return NA
    else:
        return v;

def main():
    parser = OptionParser(usage='''
        %prog [OPTIONS]
        
        .META: stdin
            1 id1
            2 id2
            3 val

        Builds a weighted incidence matrix out of a map file.
    ''')
    parser.add_option('-c', '--columns', type=str, dest='columns', default=None, help='add columns to the matrix using as their labels those specified in COLUMNS - a space separated list (in single or double quotes). Values will be NA [default: %default]', metavar='COLUMNS')
    parser.add_option('-w', '--rows', type=str, dest='rows', default=None, help='add rows to the matrix using as their labels those specified in COLUMNS - a space separated list (in single or double quotes). Values will be NA [default: %default]', metavar='ROWS')
    parser.add_option('-C', '--columns_from_file', type=str, dest='columns_from_file', default=None, help='same as -c but takes the list from file [default: %default]', metavar='COLUMNS')
    parser.add_option('-s', '--sorted', dest='sorted', action='store_true', default=False, help='assume the input block sorted on the first column (memory efficient) [default: %default]')
    parser.add_option('-S', '--sorted_column_out', dest='sorted_column_out', action='store_false', default=True, help='columns are reported in lexicographic order [default: %default]')
    parser.add_option('-k', '--kill', dest='kill', action='store_true', default=False, help='kill columns not reported in -c or -C [default: %default]')
    parser.add_option('-i', '--assume_sorted', dest='assume_sorted', action='store_true', default=False, help='do not check for lexicographic sorting of first column [default: %default]')
    parser.add_option('-t', '--transpose', dest='transpose', action='store_true', default=False, help='transpose the output matrix [default: %default]')
    parser.add_option('-e', '--missing_values', type=str, dest='missing_values', default="NA", help='specify the string used in place of missing values [default: %default]')
    parser.add_option('-r', '--row_id', type=str, dest='row_id', default=">ROW_ID", help='the label used for the first column in the header [default: %default]')
    options, args = parser.parse_args()

    if len(args) != 0:
        exit('Unexpected argument number.')
    

    if options.sorted and not (options.columns or options.columns_from_file):
        exit("In -s mode you must specify all the columns with -c.")

    if options.columns and options.columns_from_file:
        exit("-c an -C are not compatible")

    global NA
    NA=options.missing_values
    global ROW_ID
    ROW_ID=options.row_id
    
    columns = []
    if options.columns:
        columns=options.columns.split()
    elif options.columns_from_file:
        with file(options.columns_from_file, 'r') as fd:
            for line in fd:
                tokens = line.rstrip().split('\t')
                columns+=tokens

    rows=set()
    seen_rows=set()
    if options.rows:
        rows=set(options.rows.split())


    if columns:
        matrix = SparseMatrix(OrderedSet(columns))
    else:
        matrix = SparseMatrix()    
    
    xpre=None
    header_printed=False
    for line in stdin:
        try:
            if not options.transpose:
                x,y,v = line.rstrip().split("\t")
            else:
                y,x,v = line.rstrip().split("\t")

            if options.sorted:
                if y not in columns:
                    if options.kill:
                        continue
                    else:
                        exit("Unexpected column (%s) not in (%s)." % (y,",".join(columns)))
                if xpre is not None:
                    try:
                        assert(x>=xpre) or options.assume_sorted
                    except AssertionError:
                        raise Exception("Input not sorted on column 1.")
        except ValueError:
            raise Exception("Malformed input (%s)" % line)

        ##################
        #in caso di sorted stampo al cambio di x e svuoto la matrice
        #
        if options.sorted and xpre is not None and xpre!=x:#stampo linea per linea
            if not header_printed:
                print(matrix.head(sort=options.sorted_column_out))
                header_printed = True

            for r in matrix.rows():
                print('\t'.join( [None2NA(i) for i in r] ))#non posso usare v al posto di i perche` v e` gia` usata
            matrix.clean()
        #
        ##################
        xpre=x
        matrix.set(x, y, v)
        if options.rows:
            seen_rows.add(x)

    if not header_printed:
        print(matrix.head(sort=options.sorted_column_out))

    for r in matrix.rows():
        print('\t'.join( [None2NA(v) for v in r] ))
    if options.rows:
        for r in rows - seen_rows:
            print("\t".join( [r]+ [NA for c in range(len(matrix.cols_item))]))

if __name__ == '__main__':
    main()
