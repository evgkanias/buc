#include <iostream> 
#include <sstream> 
#include <fstream> 
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <sys/time.h>

#define MAX_LINE_LENGTH 4096

using namespace std;

/* types related to the output record */
typedef vector< int > V1;
typedef vector< V1 > V2;
typedef vector< V2 > V3;

/* prototypes */
void menu();
bool roll_up(int);
bool drill_down(int);
bool slice(int, int);
void init_cube();
void commit_cube();
void display_cube_summary();

int** read_datafile(int**, char*);
void display_data(int**, int, int);
int find_cardinality(int**,int, int);
void buc(int**, int, int, string);
int** partition(int**, int, int, int*, int);
void write_ancestors_one(string, int*);
void init_output();
void calculate_name(string);
void subtree_name(string, vector< string >);
int output_indexof(string);
void display_output_summary();
void commit_output();
int stoi(const string&);
string itos(int);
void elapsed_time(timeval*);

/* global vars */
int cur_cube; // This is a pointer of the current cube and output name

int minsup;
int *n_dimensions;
int *n_tuples;
V3 output;
vector< string > output_names;
timeval *prev_time;

int main(int argc, char **argv)
{
	if (argc != 3) {
		cout << "1866-1877 <datafile> <minsup>" << endl;
		exit(1);
	}

	cout << "Initializing";

	// mark start time
	prev_time = new timeval();
	gettimeofday(prev_time, NULL);

	int **data = NULL;

	/* allocate memory and init global vars */
	n_tuples = new int(0);
	n_dimensions = new int(0);
	minsup = atoi(argv[2]);

	/* read data file and get the base cuboid */
	data = read_datafile(data, argv[1]);

	/* setup output record */
	init_output();

	/* start the BUC */
	buc(data, 0, *n_tuples, "");
        init_cube();

	/* display current cube and save to disk */
        cout << endl;
        display_cube_summary();
	commit_cube();

        /* start menu */
        menu();

	/* clean up */
	delete n_dimensions;
	delete n_tuples;

	return 0;
}

void menu() {
    int choice = 0;
    int max_dim = output_names[cur_cube].size()/2;
    int dim = max_dim;

    /* prints the choices */
    do {
        cout << endl;
        cout << "Choose method:" << endl;
        cout << "1. Roll-up" << endl;
        cout << "2. Drill-down" << endl;
        cout << "3. Slice" << endl;
        cout << "4. Print the current cube" << endl;
        cout << "5. Display current cubes sammary" << endl;
        cout << "0. Cancel" << endl;
        cout << endl;
        cout << "Enter your choice > ";
        cin >> choice;
        cout << endl;

        int dimension = -1; // set an invalid dimension number
        int temp;
        int cur_dim_size = output_names[cur_cube].size()/2;
        int * cur_dim = new int[cur_dim_size];
        int i = 0;
        istringstream iss(output_names[cur_cube]);

        /* fill the current dimentions array */
        while (iss >> temp) {
            cur_dim[i] = temp + 1;
            i++;
        }

        /* print the suitable message and take the users choice */
        if (choice == 1 || choice == 2 || choice == 3) {
            if (choice == 2 && dim == max_dim) {
                cout << "You can't apply this method yet." << endl;
                cout << "Please try an other method." << endl;
                cout << endl;
                continue;
            }
            cout << "Type a dimention from the list that you want to apply the method. ";
            cout << "{ ";
            if (choice == 1 || choice == 3) {
                for (int i = 0; i < cur_dim_size; i++) {
                    cout << cur_dim[i] << " ";
                }
            } else {
                for (int i = 0; i < *n_dimensions; i++) {
                    int val = i + 1;
                    bool not_ok = false;
                    for (int j = 0; j < cur_dim_size; j++) {
                        if (val == cur_dim[j])
                            not_ok = true;
                    }
                    if (!not_ok)
                        cout << val << " ";
                }
            }
            cout << "} > ";

            cin >> dimension;
        }

        /* execute the method that the user choose */
        switch (choice) {
            case 1:
                if (dim > 1) {
                    dim--;

                    if (roll_up(dimension)) {
                        cout << "Roll-up method completed successfully!" << endl;
                        cout << "Check the output file for details." << endl;
                    }
                } else {
                    cerr << "Fail to roll-up the cube. Not enough dimentions." << endl;
                    cerr << "You have to do a drill-down first." << endl;
                }
                break;
            case 2:
                if (dim < max_dim) {
                    dim++;

                    if (drill_down(dimension)) {
                        cout << "Drill-down method completed successfully!" << endl;
                        cout << "Check the output file for details." << endl;
                    }
                } else {
                    cerr << "Fail to drill-down the cube. Too many dimentions." << endl;
                    cerr << "You have to do a roll-up first." << endl;
                }
                break;
            case 3:
                int numb;
                cout << "type the number (1-12) that you want to apply the method > ";
                cin >> numb;
                if ((dimension > 0 && dimension <= 6) && (numb > 0 && numb <= 12)) {
                    if (slice(dimension, numb)) {
                        cout << "Slice method completed successfully!" << endl;
                        cout << "Check the output file for details." << endl;
                    }
                } else {
                    cerr << "Wrong input dimension or number. Fail to slice the cube." << endl;
                    cerr << "Please try again." << endl;
                }
                break;
            case 4:
                commit_cube();
                cout << "Cube has been printed to the output file." << endl;
                cout << "Check the output file." << endl;
                break;
            case 5:
                cout << endl;
                display_cube_summary();
                cout << endl;
            case 0:
                cout << "Thank you. Bye!" << endl;
                break;
            default:
                cerr << "Wrong choice. Please try again." << endl;
        }
        cout << endl;
    } while ( choice != 0 );
}

bool roll_up(int dim) {
    dim--; // decrease the input dimension to make it compatible with our data

    /* start the roll-up algorithm */
    if (output[cur_cube].size() > 0) {
        int temp, nd_size = output_names[cur_cube].size()/2 - 1;
        int * new_dim = new int[nd_size];
        bool ok = false;
        istringstream iss(output_names[cur_cube]);
        int i = 0;

        /* fill the new dimmensions array */
        while (iss >> temp) {
            /* we skip the requested dimension */
            if (dim == temp)
                ok = true;
            else {
                new_dim[i] = temp;
                i++;
            }
        }

        if (!ok) {
            cout << "This dimension does not exist in our current cube." << endl;
            cout << "Try an other dimension." << endl;
            return false;
        }

        /* searching for the new cuboid */
        for (i = 0; i < (int) pow(2.0, (double) *n_dimensions); i++) {
            int cur_size = 0;
            istringstream iss(output_names[i]);

            /* check every cuboid (by its name) if all of its elements are
             * matched with the elements in our current dimmension array */
            while (iss >> temp) {
                ok = false;

                for (int j = 0; j < nd_size; j++) {
                    if (temp == new_dim[j])
                        ok = true;
                }
                if (!ok)
                    break;

                cur_size++;
            }

            /* takes the result of the previous step and if it's okay, checks
             * if the current dimensions array size is equal to the cuboid's
             * name size */
            if (cur_size == nd_size && ok) {
                /* sets the cube pointer to the new result */
                cur_cube = i;
                break;
            }
        }
    }

    commit_cube();
    return true;
}

bool drill_down(int dim) {
    dim--; // decrease the input dimension to make it compatible with our data

    /* start the drill-down algorithm */
    if (output[cur_cube].size() > 0) {
        int temp, nd_size = output_names[cur_cube].size()/2 + 1;
        int * new_dim = new int[nd_size];
        bool ok = true;
        istringstream iss(output_names[cur_cube]);
        int i = 0;

        /* fill the new dimmensions array */
        while (iss >> temp) {
            if (dim == temp)
                ok = false;
            else {
                new_dim[i] = temp;
                i++;
            }
        }
        /* we add the requested dimension to the array */
        new_dim[i] = dim;

        if (!ok) {
            cout << "This dimension does exist in our current cube." << endl;
            cout << "Try an other dimension." << endl;
            return false;
        }

        /* searching for the new cuboid */
        for (i = 0; i < (int) pow(2.0, (double) *n_dimensions); i++) {
            int cur_size = 0;
            istringstream iss(output_names[i]);

            /* check every cuboid (by its name) if all of its elements are
             * matched with the elements in our current dimmension array */
            while (iss >> temp) {
                ok = false;

                for (int j = 0; j < nd_size; j++) {
                    if (temp == new_dim[j])
                        ok = true;
                }
                if (!ok)
                    break;

                cur_size++;
            }

            /* takes the result of the previous step and if it's okay, checks
             * if the current dimensions array size is equal to the cuboid's
             * name size */
            if (cur_size == nd_size && ok) {
                /* sets the cube pointer to the new result */
                cur_cube = i;
                break;
            }
        }
    }

    commit_cube();
    return true;
}

bool slice(int dim, int numb) {
    /* open output file */
    ofstream outfile("output", ios::out);
    if (!outfile) {
            cout << "Error opening output file" << endl;
            exit(1);
    }
    dim--; // decrease the input dimension to make it compatible with our data

    /* start the slice algorithm */
    if (output[cur_cube].size() > 0) {
        bool ok = false;
        int p = 0, idx;
        int temp;
        istringstream iss1(output_names[cur_cube]), iss2(output_names[cur_cube]);

        /* gets the position of our dimension in the cuboid */
        while (iss1 >> temp) {
            if (dim == temp) {
                ok = true;
                idx = p;
            }
            p++;
        }

        if (!ok) {
            cout << "This dimension does not exist in our current cube." << endl;
            cout << "Try an other dimension." << endl;
            return false;
        }

        /* prints the name of the current cuboid */
        outfile << "Slice of dimension ";
        outfile << (char) (dim + 97);
        outfile << " corresponding to ";
        outfile << numb << ":" << endl;

        /* prints the result of slice method */
        for (unsigned int j = 0; j < output[cur_cube].size(); j++) {
            
            if (output[cur_cube][j][idx] != numb)
                continue;

            for (unsigned int k = 0; k < output[cur_cube][j].size(); k++) {
                outfile << output[cur_cube][j][k] << " ";
            }

            outfile << endl;
        }
        outfile << endl;
    }

    outfile.close();
    return true;
}

/* sets the current cuboid pointer
 * to the cuboid with the larger number of dimensions */
void init_cube() {
    int max = 0, counter;

    for (int i = 0; i < (int) pow(2.0, (double) *n_dimensions); i++) {
        counter = 0;
        if (output[i].size() > 0) {

            istringstream iss(output_names[i]);
            int temp;

            while (iss >> temp) {
                counter++;
            }

            if (counter > max) {
                cur_cube = i;
                max = counter;
            }
        }
    }
}

void commit_cube()
{
	/* open output file */
	ofstream outfile("output", ios::out);
	if (!outfile) {
            cout << "Error opening output file" << endl;
            exit(1);
	}

	/* write */
        if (output[cur_cube].size() > 0) {
            int temp;
            istringstream iss(output_names[cur_cube]);
            while (iss >> temp) {
                outfile << (char) (temp + 97);
            }
            outfile << endl;
            for (unsigned int j = 0; j < output[cur_cube].size(); j++) {
                for (unsigned int k = 0; k < output[cur_cube][j].size(); k++) {
                    outfile << output[cur_cube][j][k] << " ";
                }
                outfile << endl;
            }
            outfile << endl;
        }

	outfile.close();
}

void display_cube_summary()
{
    if (output[cur_cube].size() > 0) {
        int temp;
        istringstream iss(output_names[cur_cube]);
        while (iss >> temp) {
            cout << (char) (temp + 97);
        }
        cout << ":" << output[cur_cube].size() << " " << endl;
    }
}

void buc(int **data, int dim, int data_dim_size, string path)
{
	/* optimization.  if data size is 1, no need to compute.  set all
	 * ancestor's count to 1 */
	if (data_dim_size == 1) {
		write_ancestors_one(path, data[0]);
		return;
	}

	int **frequency;
	frequency = new int*[*n_dimensions];

	/* write the count to the output record */
	if ((data_dim_size > 0) & (path.length() > 0)) {
		vector< int > v;
		int temp;
		istringstream is_path(path);

		/* while there are still more tokens in the path */
		while (is_path >> temp) {
			v.push_back(data[0][temp]);
		}
		v.push_back(data_dim_size);

		int output_index = output_indexof(path);
		output[output_index].push_back(v);
	}

	for (int d = dim; d < *n_dimensions; d++) {
		/* attach the current cuboid name onto the path */
		if (path.length() == 0) { path = itos(d) + " "; }
		else { path += itos(d) + " "; }

		/* find cardinality of data on dimension d */
		int C = find_cardinality(data, d, data_dim_size);

		/* partition the data on dimension d.  also calculate the
		 * frequency of each number (dataCount in the paper) */
		frequency[d] = new int[C];
		data = partition(data, d, C, frequency[d], data_dim_size);

		/* for each partition */
		for (int i = 0; i < C; i++) {
			int c = frequency[d][i];

			/* the BUC ends here if minsup is not satisfied */
			if (c >= minsup) {
				/* construct new data set that is a subset of the
				 * original data.  this is a partition. */
				int** sub_data = new int*[c];

				/* figure out where in the original data does the
				 * particular partition start */
				int c_start = 0;
				for (int j = 0; j < i; j++) {
					c_start += frequency[d][j];
				}

				/* copy the values of the data into the new subset.
				 * note this is just a copy of the pointers.  they're
				 * still using the same actual data */
				for (int j = 0; j < c; j++) {
					sub_data[j] = data[c_start + j];
				}

				/* recursively call buc with partition */
				buc(sub_data, d + 1, c, path); 
			}
		}

		/* remove the last cuboid name in the path */
		if (path.length() == 2) { path = ""; }
		else { 
			if (path[path.length() - 3] == ' ') {
				path = path.substr(0, path.length() - 2); 
			}
			else {
				path = path.substr(0, path.length() - 3); 
			}
		}
	}


	/* clean up */
	for (int i = dim; i < *n_dimensions; i++) { delete[] frequency[i]; }
	delete[] frequency;
	if (path == "") {
		/* this completely removes the input data from memory.  it is
		 * called when the entire BUC algorith ends */
		for (int i = 0; i < *n_tuples; i++) { delete[] data[i]; }
	}
	delete[] data;
}

int** partition(int **data, int dimension, int cardinality, int
		*frequency, int data_dim_size)
{
	int *counting_sort_freq = new int[cardinality];

	/* clear frequency */
	for (int i = 0; i < cardinality; i++) { frequency[i] = 0; }

	/* calculate frequency */
	for (int i = 0; i < data_dim_size; i++) {
		frequency[data[i][dimension]]++; 
	}

	/* make copy of frequency array for counting sort */
	counting_sort_freq = new int[cardinality];
	for (int i = 0; i < cardinality; i++) { 
		counting_sort_freq[i] = frequency[i];
	}

	/* add previous frequency count to current.  i.e.
	 * counting_sort_freq[i] will mean there are counting_sort_freq[i]
	 * numbers < i */
	for (int i = 1; i < cardinality; i++) { 
		counting_sort_freq[i] += counting_sort_freq[i - 1];
	}

	/* allocate space for sorted results */
	int **sorted_data = new int*[data_dim_size];

	/* counting sort */
	for (int i = 0; i < data_dim_size; i++) {
		sorted_data[counting_sort_freq[data[i][dimension]] - 1] =
			data[i];
		counting_sort_freq[data[i][dimension]]--;
	}

	/* clean up */
	delete[] data;
	delete[] counting_sort_freq;

	return sorted_data;
}

/* returns the cardinality of data on a specific dimension */
int find_cardinality(int **data, int dimension, int data_dim_size) 
{
	/* this shouldn't happen */
	if (data_dim_size == 0) { return 0; }

	int max = data[0][dimension];
	for (int i = 1; i < data_dim_size; i++) {
		if (data[i][dimension] > max) { max = data[i][dimension]; } 
	}

	/* max + 1 because i'm assuming the range is 0..max.  this is
	 * unnecessary in this assignment because the generated data has
	 * range 1..max. */
	return (max + 1);
}

int** read_datafile(int **data, char *filename) 
{
	FILE *f;
	int temp, i, j;
	char line[MAX_LINE_LENGTH], *word;

	// open data file
	if ((f = fopen(filename, "r")) == NULL) {
		printf("Error: cannot open file %s.\n", filename);
		exit(-1);
	}

        *n_tuples = 100000;
        *n_dimensions = 6;

	/* create appropriate sized array to hold all data */
	data = new int*[*n_tuples];
	for (i = 0; i < *n_tuples; i++) { data[i] = new int[*n_dimensions]; }

	elapsed_time(prev_time);
	cout << "Reading input file";

	/* read in all data */
	i = 0;
	j = 0;

	for (i = 0; i < *n_tuples; i++) {

		// read a line
		fgets(line, MAX_LINE_LENGTH, f);

		word = strtok(line, ";");

		for (j = 0; j < *n_dimensions; j++) {
			temp = atoi(word);
			data[i][j] = temp;
			word = strtok(NULL, ";");
		}
	}

	elapsed_time(prev_time);

	/* done with all file reading.  close data file */
	fclose(f);

	return data;
}

void display_data(int **data, int dimension, int n_tuples)
{
	/* display data */
	for (int i = 0; i < n_tuples; i++) {
		for (int j = 0; j < dimension; j++) {
			cout << data[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

/* initialize output record */
void init_output()
{
	/* initialize all the data cuboids' names.  there should be exactly
	 * 2^n_dimension of these, including the "all" cuboid. */
	output_names.push_back("ALL");
	for (int i = 0; i < *n_dimensions; i++) {
		calculate_name(itos(i) + " ");
	}

	/* resize the number of data cuboids to 2^n_dimensions */
	output.resize(output_names.size());
}

/* recursively compute the parse tree */
void calculate_name(string cd)
{
	output_names.push_back(cd);

	/* find out the last cuboid in the cd string */
	int last;
	if (cd.length() == 2) {
		last = stoi(cd.substr(cd.length() - 2, 1));
	}
	else {
		if (cd[cd.length() - 3] == ' ') {
			last = stoi(cd.substr(cd.length() - 2, 1));
		}
		else {
			last = stoi(cd.substr(cd.length() - 3, 2)); 
		}
	}

	/* see if we're on the last possible dimension */
	if (last + 1 == *n_dimensions)
		return;
	else {
		for (int i = last + 1; i < *n_dimensions; i++) {
			string next_cd = cd + itos(i) + " ";
			calculate_name(next_cd);
		}
	}
}

/* returns index of a cuboid in the output record based on its name */
int output_indexof(string path)
{
	for (int i = 0; i < (int) pow(2.0, (double) *n_dimensions); i++) {
		if (output_names[i] == path) {
			return i;
		}
	}

	/* should never get here */
	return 0;
}

/* optimization. if current partition has size 1, this sets all
 * ancestors or more aggregated cuboids to 1 */
void write_ancestors_one(string path, int *tuple)
{
	vector< int > v;
	int temp;
	istringstream is_path(path);

	/* while there are still more tokens in the path */
	while (is_path >> temp) {
		v.push_back(tuple[temp]);
	}
	/* puts the count = 1 */
	v.push_back(1);

	int output_index = output_indexof(path);
	output[output_index].push_back(v);

	/* find out the last cuboid in the cd string */
	int last;
	if (path.length() == 2) {
		last = stoi(path.substr(path.length() - 2, 1));
	}
	else {
		if (path[path.length() - 3] == ' ') {
			last = stoi(path.substr(path.length() - 2, 1));
		}
		else {
			last = stoi(path.substr(path.length() - 3, 2)); 
		}
	}

	/* see if we're on the last possible dimension */
	if (last + 1 == *n_dimensions)
		return;
	else {
		for (int i = last + 1; i < *n_dimensions; i++) {
			string next_path = path + itos(i) + " ";
			write_ancestors_one(next_path, tuple);
		}
	}
}

void display_output_summary()
{
	for (int i = 0; i < (int) pow(2.0, (double) *n_dimensions); i++) {
		if (output[i].size() > 0) {
			int temp;
			istringstream iss(output_names[i]);
			while (iss >> temp) {
				cout << (char) (temp + 97);
			}
			cout << ":" << output[i].size() << " " << endl;
		}
	}
}

void commit_output()
{
	/* open output file */
	ofstream outfile("output", ios::out);
	if (!outfile) {
		cout << "Error opening output file" << endl;	
		exit(1);
	}

	/* write */
	for (int i = 0; i < (int) pow(2.0, (double) *n_dimensions); i++) {
		if (output[i].size() > 0) {
			int temp;
			istringstream iss(output_names[i]);
			while (iss >> temp) {
				outfile << (char) (temp + 97);
			}
			outfile << endl;
			for (unsigned int j = 0; j < output[i].size(); j++) {
				for (unsigned int k = 0; k < output[i][j].size(); k++) {
					outfile << output[i][j][k] << " ";
				}
				outfile << endl;
			}
			outfile << endl;
		}
	}

	outfile.close();
}

/* converts a string to an integer */
int stoi(const string &s)
{
	int result;
	istringstream(s) >> result;
	return result;
}

/* converts an integer to a string */
string itos(int n)
{
	ostringstream o;
	o << n;
	return o.str();
}

void elapsed_time(timeval *st) 
{
	timeval *et = new timeval();
	gettimeofday(et, NULL);
	long time = 1000 * (et->tv_sec - st->tv_sec) + (et->tv_usec -
			st->tv_usec)/1000;
	cout << " ... used time: " << time << " ms." << endl;
	memcpy(st, et, sizeof(timeval));
}

