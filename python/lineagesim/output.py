# -*- coding: utf-8 -*-

"""Functions for outputting data and graph structure information"""


import csv
import json


def graph_write_json(graph, filename, **kwargs):
    """Write the graph to JSON file"""
    d = {}
    d['attributes'] = {a:graph[a] for a in graph.attributes()}
# Unfortunately, doing it this way raises some errors with attributes like
# numpy arrays
#    d['vertices'] = [{'index': v.index,
#                      'attributes': v.attributes()} for v in graph.vs]
    d['vertices'] = [{'index': v.index,
                      'attributes': {'name': v['name'],
                                     'parent': v['parent'],
                                     'first_seen': v['first_seen'],
                                     'last_seen': v['last_seen'],
                                     'lineage_last_seen': v['lineage_last_seen'],
                                     'depth': v['depth'],
                                     'abundance': int(v['abundance']),
                                     'abundances': v['abundances'],
                                     'lineage_abundance': int(v['lineage_abundance']),
                                     'frequency': v['frequency'],
                                     'max_frequency': v['max_frequency'],
                                     'lineage_frequency': v['lineage_frequency'],
                                     'fixation_time': v['fixation_time'],
                                     'fitness': v['fitness'],
                                     'fitness_effects': v['fitness_effects'],
                                     'fitness_diff': v['fitness_diff']}} for v in graph.vs]

    d['edges'] = [{'index': e.index,
                   'source': e.source,
                   'target': e.target,
                   'attributes': e.attributes()} for e in graph.es]

    with open(filename, 'w') as outfile:
        json.dump(d, outfile, **kwargs)


def graph_write_csv(graph, filename, **kwargs):
    """Write the graph to CSV file"""

    with open(filename, 'w') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['Genotype',
                                                     'Parent',
                                                     'FirstSeen',
                                                     'LastSeen',
                                                     'LineageLastSeen',
                                                     'Depth',
                                                     'Fitness',
                                                     'FitnessDiff',
                                                     'NumMutations',
                                                     'Abundance',
                                                     'LineageAbundance',
                                                     'Frequency',
                                                     'MaxFrequency',
                                                     'LineageFrequency',
                                                     'FixationTime'], **kwargs)
        writer.writeheader()

        for g in graph.vs.select(lambda v: v['first_seen'] is not None):
            writer.writerow({'Genotype': g['name'],
                             'Parent': g['parent'],
                             'FirstSeen': g['first_seen'],
                             'LastSeen': g['last_seen'],
                             'LineageLastSeen': g['lineage_last_seen'],
                             'Depth': g['depth'],
                             'Fitness': g['fitness'],
                             'FitnessDiff': g['fitness_diff'],
                             'NumMutations': len(g['fitness_effects']),
                             'Abundance': g['abundance'],
                             'LineageAbundance': g['lineage_abundance'],
                             'Frequency': g['frequency'],
                             'MaxFrequency': g['max_frequency'],
                             'LineageFrequency': g['lineage_frequency'],
                             'FixationTime': g['fixation_time']})

