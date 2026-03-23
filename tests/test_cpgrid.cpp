#include <config.h>

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <dune/common/unused.hh>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <opm/grid/cpgrid/dgfparser.hh>

#if HAVE_OPM_COMMON
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#endif

#define DISABLE_DEPRECATED_METHOD_CHECK 1
using Dune::referenceElement; //grid check assume usage of Dune::Geometry
#include <dune/grid/test/gridcheck.hh>

// Re-enable warnings.
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include <iostream>

template <class GridView>
void testGridIteration( const GridView& gridView, const int nElem )
{
    typedef typename GridView::template Codim<0>::Iterator ElemIterator;
    typedef typename GridView::IntersectionIterator IsIt;
    typedef typename GridView::template Codim<0>::Geometry Geometry;

    int numElem = 0;
    ElemIterator elemIt = gridView.template begin<0>();
    ElemIterator elemEndIt = gridView.template end<0>();
    for (; elemIt != elemEndIt; ++elemIt) {
        const Geometry& elemGeom = elemIt->geometry();
        if (std::abs(elemGeom.volume() - 1.0) > 1e-8)
            std::cout << "element's " << numElem << " volume is wrong:" << elemGeom.volume() << "\n";

        typename Geometry::LocalCoordinate local( 0.5 );
        typename Geometry::GlobalCoordinate global = elemGeom.global( local );
        typename Geometry::GlobalCoordinate center = elemGeom.center();
        {
            std::cout << "center = " << center << " global( localCenter ) = " << global << std::endl;
        }

        int numIs = 0;
        IsIt isIt = gridView.ibegin(*elemIt);
        IsIt isEndIt = gridView.iend(*elemIt);
        for (; isIt != isEndIt; ++isIt, ++numIs)
        {
            const auto& intersection = *isIt;
            const auto& isGeom = intersection.geometry();
            if (std::abs(isGeom.volume() - 1.0) > 1e-8)
                std::cout << "volume of intersection " << numIs << " of element " << numElem
                          << " volume is wrong: " << isGeom.volume() << "\n";

            if (intersection.neighbor())
            {
                if (numIs != intersection.indexInInside())
                    std::cout << "num iit = " << numIs
                              << " indexInInside " << intersection.indexInInside() << std::endl;

                if (std::abs(intersection.outside().geometry().volume() - 1.0) > 1e-8)
                    std::cout << "outside element volume of intersection " << numIs
                              << " of element " << numElem
                              << " volume is wrong: "
                              << intersection.outside().geometry().volume() << std::endl;

                if (std::abs(intersection.inside().geometry().volume() - 1.0) > 1e-8)
                    std::cout << "inside element volume of intersection " << numIs
                              << " of element " << numElem
                              << " volume is wrong: "
                              << intersection.inside().geometry().volume() << std::endl;
            }
        }

        if (numIs != 2 * GridView::dimension)
            std::cout << "number of intersections is wrong for element " << numElem << "\n";

        ++numElem;
    }

    if (numElem != nElem)
        std::cout << "number of elements is wrong: " << numElem
                  << ", expected " << nElem << std::endl;
}


template <class Grid>
void testGrid(Grid& grid, const std::string& name, const size_t nElem, const size_t nVertices)
{
    typedef typename Grid::LeafGridView GridView;

    std::cout << name << std::endl;

    testGridIteration( grid.leafGridView(), nElem );

    std::cout << "create vertex mapper\n";
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper(grid.leafGridView(),
                                                               Dune::mcmgVertexLayout());
    std::cout << "VertexMapper.size(): " << mapper.size() << "\n";
    if (static_cast<size_t>(mapper.size()) != nVertices)
        std::cout << "Wrong size of vertex mapper. Expected " << nVertices << "!" << std::endl;

    if (true || grid.geomTypes(0)[0].isCube())
    {
        std::cout << "create vtkWriter\n";
        typedef Dune::VTKWriter<GridView> VtkWriter;
        VtkWriter vtkWriter(grid.leafGridView());

        std::cout << "create cellData\n";
        int numElems = grid.size(0);
        std::vector<double> tmpData(numElems, 0.0);

        std::cout << "add cellData\n";
        vtkWriter.addCellData(tmpData, name);

        std::cout << "write data\n";
        vtkWriter.write(name, Dune::VTK::ascii);
    }
}


    // In test_cpgrid.cpp or a utility header:
std::pair<std::string, std::string> splitDeck(const std::string& combined,
                                               const std::string& delimiter = "-- GRID2")
{
    auto pos = combined.find(delimiter);
    if (pos == std::string::npos)
        throw std::runtime_error("Delimiter '" + delimiter + "' not found in deck string");
    return { combined.substr(0, pos),
             combined.substr(pos + delimiter.size()) };
}



int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    typedef Dune::CpGrid Grid;

#if HAVE_OPM_COMMON



    const char* deckStringcomposite =
        "-- GRID1\n"
        "RUNSPEC\n"
        "METRIC\n"
        "DIMENS\n"
        "2 2 2 /\n"
        "GRID\n"
        "DXV\n"
        "2*1 /\n"
        "DYV\n"
        "2*1 /\n"
        "DZ\n"
        "8*1 /\n"
        "TOPS\n"
        "8*100.0 /\n"
        "PORO\n"
        "8*0.8 /\n"
        "PERMX\n"
        "8*6 /\n"
        "-- GRID2\n"           //  delimiter gri2
        "RUNSPEC\n"
        "METRIC\n"
        "DIMENS\n"
        "3 3 3 /\n"
        "GRID\n"
        "DXV\n"
        "3*1 /\n"
        "DYV\n"
        "3*1 /\n"
        "DZ\n"
        "27*1 /\n"
        "TOPS\n"
        "27*100.0 /\n"
        "PORO\n"
        "27*0.3 /\n"
        "PERMX\n"
        "27*10 /\n";


    const char *deckString =
        "RUNSPEC\n"
        "METRIC\n"
        "DIMENS\n"
        "2 2 2 /\n"
        "GRID\n"
        "DXV\n"
        "2*1 /\n"
        "DYV\n"
        "2*1 /\n"
        "DZ\n"
        "8*1 /\n"
        "TOPS\n"
        "8*100.0 /\n"
        "PORO\n"
        "8*0.2 /\n"
        "PERMX\n"
        "8*5 /\n";

    const char *deckString2 =
        "RUNSPEC\n"
        "METRIC\n"
        "DIMENS\n"
        "3 3 3 /\n"
        "GRID\n"
        "DXV\n"
        "3*1 /\n"
        "DYV\n"
        "3*1 /\n"
        "DZ\n"
        "27*1 /\n"
        "TOPS\n"
        "27*100.0 /\n"
        "PORO\n"
        "27*0.3 /\n"
        "PERMX\n"
        "27*10 /\n";

    Opm::Parser parser;
    // const auto deck  = parser.parseString(deckString);
    // const auto deck2 = parser.parseString(deckString2);




// Usage:
auto [deckStr1, deckStr2] = splitDeck(deckStringcomposite);
const auto deck  = parser.parseString(deckStr1);
const auto deck2 = parser.parseString(deckStr2);



    Opm::EclipseGrid ecl_grid(deck);
    Opm::EclipseGrid ecl_grid2(deck2);
    Opm::TableManager tables(deck);
    Opm::TableManager tables2(deck2);


    Grid grid;
    Grid grid2;

    grid.processEclipseFormat(&ecl_grid,  nullptr, false, false, false);
    grid2.processEclipseFormat(&ecl_grid2, nullptr, false, false, false);

    testGrid(grid,  "CpGrid_ecl",  8,  27);
    testGrid(grid2, "CpGrid_ecl2", 27, 64);

    // const int grid1_nc = grid.numCells();  //  8

    // grid.extendGrid(grid2, {0.0, 0.0, 2.5}, false);
    // testGrid(grid, "CpGrid_eclcombined", 35, 91);

    // Opm::FieldPropsManager fp1(deck,  Opm::Phases{true, true, true}, ecl_grid,  tables);
    // Opm::FieldPropsManager fp2(deck2, Opm::Phases{true, true, true}, ecl_grid2, tables2);

    // const auto& poro1 = fp1.get_double("PORO");  
    // const auto& poro2 = fp2.get_double("PORO");  
    // const int nCells = grid.numCells();  // 35 after extension
    // std::vector<double> merged_poro(nCells);

    // for (int c = 0; c < nCells; ++c) {
    //     int cartIdx = grid.globalCell()[c];
    //     if (grid.cellGridOrigin(c) == 0) {
    //         // Cell from grid1: cartIdx is 0-based into ecl_grid (0..7)
    //         merged_poro[c] = poro1[cartIdx];
    //     } else {
    //         // Cell from grid2: cartIdx has been offset by grid1_nc, subtract it back
    //         merged_poro[c] = poro2[cartIdx - grid1_nc];
    //     }
    // }



        
const int grid1_nc = grid.numCells();
grid.extendGrid(grid2, {0.0, 0.0, 2.5}, false);

Opm::FieldPropsManager fp1(deck,  Opm::Phases{true,true,true}, ecl_grid,  tables);
Opm::FieldPropsManager fp2(deck2, Opm::Phases{true,true,true}, ecl_grid2, tables2);

fp1.extendGrid(fp2);  // fp1 now covers all 35 cells

const auto& poro = fp1.get_double("PORO");  // size = 35
const auto& permx = fp1.get_double("PERMX");  // size = 35

for (int c = 0; c < grid.numCells(); ++c) {
    int cartIdx = grid.globalCell()[c];
    std::cout << "cell " << c
              << "  origin = "   << grid.cellGridOrigin(c)
              << "  centroid = " << grid.cellCentroid(c)
              << "  poro = "     << poro[cartIdx]
              << "  permx = "     << permx[cartIdx]
              << "\n";
}

std::cout << "create vtkWriter\n";
typedef typename Grid::LeafGridView GridView;
typedef Dune::VTKWriter<GridView> VtkWriter;
VtkWriter vtkWriter(grid.leafGridView());

std::cout << "create cellData\n";
int numElems = grid.size(0);
std::vector<double> tmpData(numElems, 0.0);

std::cout << "add cellData\n";
vtkWriter.addCellData(poro, "poro");
vtkWriter.addCellData(permx, "permx");

std::cout << "write data\n";
vtkWriter.write("outputfile", Dune::VTK::ascii);


    // std::cout << "\n--- Per-cell property summary ---\n";
    // for (int c = 0; c < nCells; ++c) {
    //     std::cout << "cell " << c
    //               << "  cartIdx = "  << grid.globalCell()[c]
    //               << "  origin = "   << grid.cellGridOrigin(c)
    //               << "  centroid = " << grid.cellCentroid(c)
    //               << "  poro = "     << merged_poro[c]
    //               << "\n";
    // }

    // const auto& grid_leafView = grid.leafGridView();
    // Dune::CartesianIndexMapper<Grid> grid_cartMapper(grid);

    // std::cout << "\n--- CartesianIndexMapper check ---\n";
    // std::cout << "cart mapper index of cell 0: "
    //           << grid_cartMapper.cartesianIndex(0) << "\n";

    // for (const auto& element : Dune::elements(grid_leafView)) {
    //     const auto& centroid_entity = grid.getEclCentroid(element);
    //     const auto& centroid_index  = grid.getEclCentroid(element.index());
    //     for (int coord = 0; coord < 3; ++coord)
    //         assert(centroid_entity[coord] == centroid_index[coord]);
    // }
    // std::cout << "getEclCentroid consistency check passed.\n";

#endif

    return 0;
}