#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 12/01/2016
"""

import json
from networkx.readwrite import json_graph


def import_graph(file_name):
    with open(file_name) as f:
        data = json.load(f)
        return json_graph.adjacency_graph(data)
