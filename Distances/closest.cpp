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
    if( !( argc == 3
           || argc == 4
         )
      )
        {
            std::cout << "Wrong number of arguments!\n";
            std::cout << "USAGE: " << argv[ 0 ] << " in_fasta search_fasta (optional) num_threads\n";
            std::cout << "If search_fasta is not included, then in_fasta will be both the reference and the search.\n";
            return EXIT_FAILURE;
        }

    std::string in_f = argv[ 1 ];

    int num_threads_index = 0;

    std::vector<std::string> seqs;
    std::vector<std::string> search_seqs;
    std::vector<std::string> *search_seqs_ptr = nullptr;

    get_seqs( seqs, in_f );
	std::size_t length = seqs[ 0 ].length();

    if( argc == 3 )
        {
            num_threads_index = 2;
            search_seqs_ptr = &seqs;
        }
    else // argc = 4 here
        {
            std::string search_f = argv[ 2 ];
            num_threads_index = 3;

            get_seqs( search_seqs, search_f );
            search_seqs_ptr = &search_seqs;
        }

    omp_set_num_threads( atoi( argv[ num_threads_index ] ) );
    

    std::unordered_map<std::string, std::pair<int,std::string>> map;
    map.reserve( search_seqs_ptr->size() );

    for( const auto& seq : *search_seqs_ptr )
        {
            map.emplace( seq, std::make_pair( length, "" ) );
        }

    std::size_t index = 0;
    std::size_t inner_index = 0;

    #pragma omp parallel for private( index, inner_index ) shared( seqs, search_seqs_ptr, map )
    for( index = 0; index < search_seqs_ptr->size(); ++index )
        {
            const std::string& my_seq = (*search_seqs_ptr)[ index ];
            auto distance = map.find( my_seq );

            for( inner_index = 0; inner_index < seqs.size(); ++inner_index )
                {
                    const std::string& my_seq_inner = seqs[ inner_index ];

                    int dist = hamming( my_seq, my_seq_inner );

                    if( dist > 0
                        && distance->second.first > dist
                      )
                        {
                            distance->second.first = dist;
                            distance->second.second = my_seq_inner;
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
