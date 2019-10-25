#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <ostream>

#include <omp.h>

void get_seqs( std::vector<std::string>& seqs, std::string in_f );
int hamming( const std::string& a, const std::string& b );
    

int main( int argc, char **argv )
{
    if( argc != 3 )
        {
            std::cout << "Wrong number of arguments, must include input fasta and number of threads";
            return EXIT_FAILURE;
        }

    std::string in_f = argv[ 1 ];
    omp_set_num_threads( atoi( argv[ 2 ] ) );

    std::vector<std::string> seqs;

    get_seqs( seqs, in_f );

	std::size_t length = seqs[ 0 ].length();

    std::vector<std::vector<int>> all_pairs;
    all_pairs.reserve( seqs.size() );

    for( int x = 0; x < seqs.size(); ++x )
        {
            all_pairs.emplace_back( std::vector<int>() );
            all_pairs[ x ].reserve( seqs.size() );
        }

    std::size_t index = 0;
    std::size_t inner_index = 0;

    #pragma omp parallel for private( index, inner_index ) shared( seqs, length )
    for( index = 0; index < seqs.size(); ++index )
        {
            std::string my_seq = seqs[ index ];
            std::vector<int>& my_vec = all_pairs[ index ];

            for( inner_index = 0; inner_index < seqs.size(); ++inner_index )
                {
                    my_vec[ inner_index ] = hamming( my_seq, seqs[ inner_index ] );
                }
        }

    std::ofstream out( "distances.tsv", std::ofstream::out );

    out << "\t";

    for( int idx = 0; idx < seqs.size() - 1; ++idx )
        {
            out << seqs[ idx ] << "\t";
        }

    out << seqs[ seqs.size() - 1 ] << "\n";

    for( int index = 0; index < seqs.size(); ++index )
        {
            std::string& current = seqs[ index ]; 

            out << current << "\t";
            std::vector<int>& vec = all_pairs[ index ];

            for( int inner_index = 0; inner_index < seqs.size() - 1; ++inner_index )
                {
                    out << vec[ inner_index ] << "\t";
                }
            out << vec[ seqs.size() - 1 ] << "\n";
        }

}

void get_seqs( std::vector<std::string>& seqs, std::string in_f )
{
    std::ifstream input_file( in_f, std::ios_base::in );
    std::string line;

    while( std::getline( input_file, line ).good() )
        {
            if( line.length() > 0
                && line[ 0 ] != '>'
              )
                {
                    seqs.emplace_back( line );
                }
        }
}

int hamming( const std::string& a, const std::string& b )
{
    std::size_t index = 0;
    int dist = 0;

    for( index = 0; index < a.length(); ++index )
        {
            dist += a[ index ] != b[ index ];
        }
    return dist;
}
