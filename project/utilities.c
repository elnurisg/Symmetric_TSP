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

// her defe yaradanda, artiq qalir o arraylar??
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

void shuffleArray(int arr[], int n) {
    srand(time(NULL));
// shuffle the array with keeping segments together
// we use pair of numbers in order to indicate the segments (its starting and ending)
    for (int i = n - 1; i > 0; i -= 2) {
        int j = rand() % (i + 1);

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
