/*
 *                     GNU GENERAL PUBLIC LICENSE
 *                       Version 2, June 1991
 *
 * Copyright (C) 2016       @arlenebatada
 *
 * Copyright (C) 1989, 1991 Free Software Foundation, Inc., <http://fsf.org/>
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 * Everyone is permitted to copy and distribute verbatim copies
 * of this license document, but changing it is not allowed.
 *
 * MODIFIED: Muhammad Haseeb, 2020
 *
 */

#pragma once

#include "common.h"
#include <limits>
#include <iostream>

/* Defines */
#define ERROR_SIZE_LESS_THAN_1 -2
#define ERROR_SIZE_INCREASED -3
#define ERROR_POSITION_GREATER_THAN_HEAP_SIZE size
#define ERROR_SIZE_DECREASED -4
#define ERROR_POSITION_LESS_THAN_0 -1
#define ERROR_HEAP_FULL capacity

/* Main class */
template<class T>
class minHeap //the main min heap class
{
private:
    int size;
    int capacity;
    T* array;

    int heapify(int element_position); 
    void swap(T&, T&);
    void swap(int, int);

public:

    minHeap()
    {
        this->capacity = 0;
        this->size = 0;
        this->array = NULL;
    }

    minHeap(int capacity)
    {
        this->capacity = capacity;
        this->size = 0;
        this->array = new T[this->capacity];
    }

    ~minHeap()
    {
        this->capacity = 0;
        this->size = 0;

        if (this->array != NULL)
        {
            delete[] this->array;
            this->array = NULL;
        }
    }

    int init(int capacity);
    int reset();
    int insert(T &element); 
    int get_capacity();
    int get_size();
    T extract_min(); 
    int decrease_key(int element_position, T new_value);
    int increase_key(int element_position, T new_value);
    int heap_sort(T *output_array); 
    T show_element(int element_position);
    T getMax();
};



template<class T>
int minHeap<T>::init(int capacity)
{
    this->capacity = capacity;
    this->size = 0;
    this->array = new T[this->capacity];

    return 0;
}

template<class T>
int minHeap<T>::reset()
{
    this->size = 0;
    for (int i = 0; i < this->capacity; i++)
    {
        T fresh;
        array[i] = fresh;
    }
    return 0;
}

template<class T>
int minHeap<T>::insert(T &element)
{
    if (size == capacity)
    {
        return increase_key(0, element);
    }

    array[size++] = element;

    for (int i = (size - 1) / 2; i >= 0; i--)
    {
        heapify(i);
    }

    return 0;
}

template<class T>
int minHeap<T>::get_capacity() { return capacity; }

template<class T>
int minHeap<T>::get_size() { return size; }

template<class T>
void minHeap<T>::swap(T& t1, T& t2)
{
    T temp = t1;
    t1 = t2;
    t2 = temp;
}

template<class T>
void minHeap<T>::swap(int p1, int p2)
{
    T&& temp  = array[p1];
    array[p1] = array[p2];
    array[p2] = temp;
}

template<class T>
int minHeap<T>::heapify(int element_position)
{
    int rchild_pos = (element_position + 1) * 2;
    int lchild_pos = rchild_pos - 1;

    int smallest_element_position = element_position;

    T smallest_element = array[smallest_element_position];

    if (lchild_pos < size && array[lchild_pos] < smallest_element)
    {
        smallest_element_position = lchild_pos;
        smallest_element          = array[lchild_pos];
    }

    if (rchild_pos < size && array[rchild_pos] < smallest_element)
    {
        smallest_element_position = rchild_pos;
        smallest_element          = array[rchild_pos];
    }

    if (smallest_element_position != element_position)
    {
        swap(array[smallest_element_position], array[element_position]);
        heapify(smallest_element_position);
    }

    return 0;
}

template<class T>
T minHeap<T>::extract_min() 
{
    if (size < 1)
    {
        return ERROR_SIZE_LESS_THAN_1;
    }

    T min = array[0];
    swap(array[0], array[size - 1]);
    size--; 
    heapify(0);
    return min;
}

template<class T>
int minHeap<T>::decrease_key(int element_position, T new_value)
{
    if (new_value > array[element_position]) //if an attempt to increase the value of the element
    {
        return ERROR_SIZE_INCREASED;
    }

    array[element_position] = new_value;

    int parent_position = (element_position + 1 >> 1) - 1;

    while ((array[parent_position] > array[element_position]) && (parent_position >= 0))
    {
        swap(array[parent_position], array[element_position]);
        element_position = parent_position;
        parent_position = (parent_position - 1) >> 1;
    }
    return 0;
}

template<class T>
int minHeap<T>::increase_key(int element_position, T new_value)
{
    if (new_value < array[element_position])
    {
        return ERROR_SIZE_DECREASED;
    }

    array[element_position] = new_value;

    heapify(element_position);

    return 0;
}

template<class T>
int minHeap<T>::heap_sort(T *output_array)
{
    if (size < 1)
    {
        return ERROR_SIZE_LESS_THAN_1;
    }

    int max_loop_count = size;

    for (int i = 0; i < max_loop_count; i++)
    {
        output_array[i] = extract_min();
    }

}

template<class T>
T minHeap<T>::show_element(int element_position)
{
    if (element_position > size - 1)
    {
        return ERROR_POSITION_GREATER_THAN_HEAP_SIZE;
    }

    return array[element_position];
}

template<class T>
T minHeap<T>::getMax()
{
    if (size <= 1)
    {
        return array[0];
    }

    T maxe = array[size - 1];

    // the limit here was k > 0 instead of size/2
    for (int_t k = size-2; k >= size/2; k --)
    {
        if (array[k] > maxe)
        {
            maxe = array[k];
        }
    }

    return maxe;
}