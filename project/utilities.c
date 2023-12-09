#include <stdio.h>
#include <stdlib.h>

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