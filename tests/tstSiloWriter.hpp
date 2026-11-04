#include <Cabana_Core.hpp>
#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>
#include <ProblemManager.hpp>
#include <SiloWriter.hpp>

#include <mpi.h>
#include <silo.h>

#include <climits>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "MeshFixture.hpp"
#include "ProblemManagerFixture.hpp"
#include "tstDriver.hpp"

/*
 * The writer names its files relative to the working directory, so each
 * test runs in a fresh temporary directory that rank 0 creates and every
 * rank enters. Tests at different rank counts run concurrently.
 */
template <class T>
class SiloWriterTest : public ProblemManagerTest<T>
{
  protected:
    using writer_type = CabanaFluids::SiloWriter<2, typename T::ExecutionSpace,
                                                 typename T::MemorySpace>;
    std::filesystem::path startDir_;
    std::filesystem::path tmpDir_;

    virtual void SetUp() override
    {
        ProblemManagerTest<T>::SetUp();
        char dir[PATH_MAX] = {};
        if ( this->testPM_->mesh()->rank() == 0 )
        {
            auto pattern = std::filesystem::temp_directory_path() /
                           "cabanafluids-silo-XXXXXX";
            pattern.string().copy( dir, PATH_MAX - 1 );
            ASSERT_NE( mkdtemp( dir ), nullptr );
        }
        MPI_Bcast( dir, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD );
        tmpDir_ = dir;
        startDir_ = std::filesystem::current_path();
        std::filesystem::current_path( tmpDir_ );
    }

    virtual void TearDown() override
    {
        MPI_Barrier( MPI_COMM_WORLD );
        std::filesystem::current_path( startDir_ );
        if ( this->testPM_->mesh()->rank() == 0 )
            std::filesystem::remove_all( tmpDir_ );
        ProblemManagerTest<T>::TearDown();
    }
};

TYPED_TEST_SUITE( SiloWriterTest, MeshDeviceTypes );

TYPED_TEST( SiloWriterTest, CreatesOutputDirectories )
{
    typename TestFixture::writer_type writer( this->testPM_ );
    writer.siloWrite( "Mesh", 0, 0.0, 0.01 );
    MPI_Barrier( MPI_COMM_WORLD );

    EXPECT_TRUE( std::filesystem::exists( "data/CabanaFluids00000.silo" ) );
    EXPECT_TRUE( std::filesystem::exists(
        "data/raw/CabanaFluidsOutput0000000000.silo" ) );
}

TYPED_TEST( SiloWriterTest, RecordsCycleAndTime )
{
    const int cycle = 3;
    const double time = 1.2345;
    typename TestFixture::writer_type writer( this->testPM_ );
    writer.siloWrite( "Mesh", cycle, time, 0.01 );
    MPI_Barrier( MPI_COMM_WORLD );

    // Rank 0 is always in the first PMPIO group, so its domain is in that
    // group's file
    if ( this->testPM_->mesh()->rank() != 0 )
        return;
    DBfile* file = DBOpen( "data/raw/CabanaFluidsOutput0000000003.silo",
                           DB_UNKNOWN, DB_READ );
    ASSERT_NE( file, nullptr );
    ASSERT_EQ( DBSetDir( file, "domain_00000" ), 0 );

    DBquadmesh* mesh = DBGetQuadmesh( file, "Mesh" );
    ASSERT_NE( mesh, nullptr );
    EXPECT_EQ( mesh->cycle, cycle );
    EXPECT_FLOAT_EQ( mesh->time, static_cast<float>( time ) );
    EXPECT_DOUBLE_EQ( mesh->dtime, time );
    DBFreeQuadmesh( mesh );

    DBquadvar* quantity = DBGetQuadvar( file, "quantity" );
    ASSERT_NE( quantity, nullptr );
    EXPECT_EQ( quantity->cycle, cycle );
    EXPECT_FLOAT_EQ( quantity->time, static_cast<float>( time ) );
    EXPECT_DOUBLE_EQ( quantity->dtime, time );
    DBFreeQuadvar( quantity );

    DBClose( file );
}

TYPED_TEST( SiloWriterTest, ReportsUnwritableOutput )
{
    // A regular file named data blocks the output directories
    if ( this->testPM_->mesh()->rank() == 0 )
        std::ofstream( "data" ) << "not a directory\n";
    MPI_Barrier( MPI_COMM_WORLD );

    typename TestFixture::writer_type writer( this->testPM_ );
    EXPECT_THROW( writer.siloWrite( "Mesh", 0, 0.0, 0.01 ),
                  std::runtime_error );
}
