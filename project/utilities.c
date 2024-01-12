#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
    int last_pos = 1;

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