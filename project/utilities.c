#include "utilities.h"

void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }   

double random01() { return ((double) rand() / RAND_MAX); } // return a random value in range 0.0-1.0

double cost(int i, int j, instance *inst)
{
	return inst->cost[i*inst->nnodes + j];
}


// plot the solution in commands.txt file 
//and writes the solution to the file if writing_to_file is true
void plot_tsp_tour(instance *inst, int writing_to_file){

	if (writing_to_file == 1)
	{
		FILE *f = fopen("plot/data.dat", "w");
		if (f == NULL)
		{
			printf("Error opening file!\n");
			exit(1);
		}
		for (int i = 0; i < inst->nnodes+1; i++)
		{
			fprintf(f, "%f %f\n", inst->xcoord[inst->best_sol[i]], inst->ycoord[inst->best_sol[i]]);
		}
		fclose(f);	
	}
	

    system("gnuplot ./plot/commands.txt");

}

int verify_tour(instance *inst){
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = 0; j < inst->nnodes; j++)
		{
			if(i!=j && inst->best_sol[i]==inst->best_sol[j]) {
				printf("indexes %d and %d are the same nodes: %d\n", i,j, inst->best_sol[i]);
				return 1;
			}

		}

	}
	if (inst->best_sol[0] != inst->best_sol[inst->nnodes]){
		printf("first and the last nodes are not the same\n");
		return 1;
	}
	
	return 0;
}


int dist(int i, int j, instance *inst)
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j]; 
	int dis = sqrt(dx*dx+dy*dy) + 0.5; // nearest integer 
	return dis;
}        

int random_0_to_length_but_different_than_previous(instance *inst, int length, int previous_random_value){

	int random_value = random_0_to_length(inst, length);

	while (random_value == previous_random_value)
	{
		random_value = random_0_to_length(inst, length);
	}
	
	return random_value;
}


int random_node_with_time_seed(int length){
	srand((unsigned int)time(NULL));
	int random_position = rand() % length;

	if (random_position<0 && random_position>=length)
	{
		printf("\n Error! Position exceeds the size of array");
		return -1;
	}
	
	return random_position;
}

int random_0_to_length(instance *inst, int length){

	if (inst->random_seed !=0){
		srand(2635623+abs(inst->random_seed));
		for (size_t i = 0; i < 1000; i++) random();
	}

	int random_position = rand() % length;
	if (random_position<0 && random_position>=length)
	{
		printf("\n Error! Position exceeds the size of array");
		return -1;
	}
	
	return random_position;
}


int * copy_array(int *arr, int size){
    int *cpy = (int *)calloc(size, sizeof(int));
    for (int i = 0; i < size; i++)
    {
        cpy[i] = arr[i];
    }

    return cpy; 
}

int * add_to_array(int starting_pos, int new_element, int *arr, int size){
    int *new_arr = (int *)calloc(size, sizeof(int));

    for (int i = 0; i < size; i++)
    {
        if (i <= starting_pos)
        {
            new_arr[i] = arr[i];
        }
        else if (i == starting_pos+1)
        {
            new_arr[i] = new_element;
        }
        else
        {
            new_arr[i] = arr[i-1];
        }
        
    }

    return new_arr;

}

int * remove_from_array(int pos, int *arr, int size){
    int *new_arr = (int *)calloc(size, sizeof(int));

    arr[size-1] = -1;

    for (int i = 0; i < size-1; i++)
    {
        if (i < pos)
        {
            new_arr[i] = arr[i];
        }
        else
        {
            new_arr[i] = arr[i+1];
        }
        
    }

    new_arr[size-1] = -1;
    return new_arr;

}


void print_array(int *arr, int size){
    printf("array is: ");
    for (int i = 0; i < size; i++)
    {
        printf("index %d: %d || ", i, arr[i]);
    }
    printf("\n\n");
}


int compare(const void* num1, const void* num2) // comparing function for sorting
{  
    int a = *(int*) num1;  
    int b = *(int*) num2;  
    if(a > b)  
    {  
        return 1;  
    }  
    else if(a < b)  
    {  
        return -1;  
    }  
    return 0;  
} 


void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void shuffle_array_for_kick(int arr[], int n) {
    srand(time(NULL)); int j;
// shuffle the array with keeping segments together
// we use pair of numbers in order to indicate the segments (its starting and ending)
    for (int i = n - 1; i > 0; i -= 2) {
        j = rand() % (i + 1);

        // swap elements i and j
        swap(&arr[i], &arr[j]);
        if (j % 2 == 0) //if new position is even number then 
        {   // we should add another element of pair after
            if (i % 2 == 0)         swap(&arr[i + 1], &arr[j + 1]);
            else                    swap(&arr[i - 1], &arr[j + 1]);
    // and also we should know about old order of pairs
    // which side of the pair was the element so that we can find its pair.
            
        }
        else // if it's odd then we should add it before
        {
            if (i % 2 == 0)         swap(&arr[i + 1], &arr[j - 1]);
            else                    swap(&arr[i - 1], &arr[j - 1]);
        }
        
        
    }
}

// function to shuffle an array of any type using Fisher-Yates shuffle
void shuffle_tsp_sol(int *arr, int nnodes) {
    int j; int temp;
    for (int i = nnodes - 1; i > 0; i--) {
        j = rand() % (i + 1);

        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }

    // close the tour
    arr[nnodes] = arr[0];
}


int * combine_two_tours_from_pos(int *arr1, int *arr2, int size, int position) {

    int *result_arr = (int *) calloc(size, sizeof(int));

    // copy elements from arr1 before the point
    for (int i = 0; i < position; i++) {
        result_arr[i] = arr1[i];
    }

    // copy elements from arr2 starting from the point
    for (int i = position; i < size; i++) {
        result_arr[i] = arr2[i];
    }

    result_arr[size] = result_arr[0]; // close the tour

    return result_arr;
}

void write_cost_to_file(double cost, const char *filename, int append) {
    const char *folder = "cost_plot";
    struct stat st = {0};
    // Check if the folder exists
    if (stat(folder, &st) != 0){
        //if the folder doesn't exist, create it
        // For Linux/Unix
        if (mkdir(folder, 0777) != 0) {
            printf("Error creating folder: %s\n", folder);
            return;
        }
    }

    FILE *fp = NULL;  
    if (append == 0)
        fp = fopen(filename, "a");  // Append mode
    else if (append == 1)
        fp = fopen(filename, "w");  // Overwrite mode
    
    if (fp == NULL) {  
        printf("Error opening file.\n");
        return;
    }

    fprintf(fp, "$$%.2f\n", cost);

    fclose(fp);
}
