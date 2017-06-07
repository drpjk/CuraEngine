/*
 * Copyright (c) 2017 Ultimaker B.V.
 *
 * CuraEngine is released under the terms of the AGPLv3 or higher.
 */
#include <cassert>
#include <cmath>
#include <cstdlib>

#include <limits>
#include <algorithm>
#include <vector>

#ifndef NDEBUG
# include <iostream>
#endif // NDEBUG

#include "math.h"
#include "FermatSpiral.h"
#include "FermatSpiralMST.h"

using namespace cura;


/*
 * Starting from the given point, travel along the given path for the given distance,
 * and return the point.
 */
static ClipperLib::IntPoint get_point_on_path_after_distance(const ClipperLib::Path& path,
                                                             uint64_t start_point_index,
                                                             int64_t distance);


void FermatSpiralInfillGenerator::generateInfill(
    const Polygons& in_outline,
    Polygons& result_lines,
    const SliceMeshStorage* mesh)
{
    // construct the MST
    SpiralContourTree tree;
    tree.setPolygons(in_outline.getPaths());
    tree.constructTree();

    // TODO: connect the contours
    ClipperLib::Path full_path;
    tree.connectContours(full_path);
}
