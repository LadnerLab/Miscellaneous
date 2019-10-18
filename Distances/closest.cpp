#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <cstdlib>
#include <algorithm>

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

    std::unordered_map<std::string, std::pair<int,std::string>> map;
    map.reserve( seqs.size() );

    for( const auto& seq : seqs )
        {
            map.emplace( seq, std::make_pair(length, "") );
        }

    std::size_t index = 0;
    std::size_t inner_index = 0;

    #pragma omp parallel for private( index, inner_index ) shared( seqs, map )
    for( index = 0; index < seqs.size(); ++index )
        {
            std::string my_seq = seqs[ index ];

            for( inner_index = 0; inner_index < seqs.size(); ++inner_index )
                {
                    if( inner_index != index )
                        {
                            int dist = hamming( my_seq, seqs[ inner_index ] );

                            if( map.find( my_seq )->second.first > dist )
                                {
                                    map.find( my_seq )->second = std::make_pair( dist, seqs[ inner_index ] );
                                }
                        }
                }
        }

	std::cout << "Sequence\tDistance\tNearest Neighbor\n";
    for( const auto& it : map )
        {
			if( it.second.second.compare( "" ) )
			{
					std::cout << it.first << "\t" << it.second.first << "\t" << it.second.second << "\n";
			}
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
