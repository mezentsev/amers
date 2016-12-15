//
// Created by Vadim Mezentsev on 15.12.16.
//

#include "connectivity.h"

p4est_connectivity_t *new_connectivity(void) {
    const p4est_topidx_t num_vertices = 8;
    const p4est_topidx_t num_trees = 1;
    const p4est_topidx_t num_ett = 0;
    const p4est_topidx_t num_ctt = 0;
    const double        vertices[8 * 3] = {
            0, 0, 0,
            100, 0, 0,
            0, 100, 0,
            100, 100, 0,
            0, 0, 100,
            100, 0, 100,
            0, 100, 100,
            100, 100, 100,
    };
    const p4est_topidx_t tree_to_vertex[1 * 8] = {
            0, 1, 2, 3, 4, 5, 6, 7,
    };
    const p4est_topidx_t tree_to_tree[1 * 6] = {
            0, 0, 0, 0, 0, 0,
    };
    const int8_t        tree_to_face[1 * 6] = {
            0, 1, 2, 3, 4, 5,
    };

    return p4est_connectivity_new_copy (num_vertices, num_trees, 0, 0,
                                        vertices, tree_to_vertex,
                                        tree_to_tree, tree_to_face,
                                        NULL, &num_ett, NULL, NULL,
                                        NULL, &num_ctt, NULL, NULL);
}
