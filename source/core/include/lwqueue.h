/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
 */

#pragma once 

#include "common.hpp"
#include <semaphore.h>

using namespace std;

#define DEF_CAPACITY                      4
#define ISFULL                           -1
#define ISEMPTY                          -2

template <class T>
class lwqueue
{
private:
    T     *arr;
    int_t filled;
    int_t    cap;
    int_t   head;
    int_t   tail;
    sem_t lock;
    BOOL isSem;

public:

    lwqueue()
    {
        cap = DEF_CAPACITY;
        arr = new T[cap];
        filled = 0;
        head = 0;
        tail = -1;
        sem_init(&lock, 0, 1);
        isSem = true;
    }

    lwqueue(int_t dcap)
    {
        cap = dcap;
        arr = new T[cap];
        filled = 0;
        head = 0;
        tail = -1;
        sem_init(&lock, 0, 1);
        isSem = true;
    }

    lwqueue(BOOL sem)
    {
        cap = DEF_CAPACITY;
        arr = new T[cap];
        filled = 0;
        head = 0;
        tail = -1;

        if (sem == true)
        {
            sem_init(&lock, 0, 1);
            isSem = true;
        }
        else
        {
            isSem = false;
        }
    }

    lwqueue(int_t dcap, BOOL sem)
    {
        cap = dcap;
        arr = new T[cap];
        filled = 0;
        head = 0;
        tail = -1;
        if (sem)
        {
            sem_init(&lock, 0, 1);
            isSem = true;
        }
        else
        {
            isSem = false;
        }
    }

    lwqueue(T *ar, int_t sz, BOOL sem)
    {
        cap = sz;
        arr = ar;
        filled = sz;
        head = 0;
        tail = sz - 1;

        if (sem)
        {
            sem_init(&lock, 0, 1);
            isSem = true;
        }
        else
        {
            isSem = false;
        }
    }

    virtual ~lwqueue()
    {
        delete[] arr;
        arr = NULL;
        cap = 0;
        filled = 0;
        head = 0;
        tail = -1;

        if (isSem)
        {
            sem_destroy(&lock);
            isSem = false;
        }
    }

    inline int_t plusOne(int_t idx)
    {
        return ((idx + 1) % cap);
    }

    inline int_t minusOne(int_t idx)
    {
        return ((idx - 1) % cap);
    }

    status_t push(T elmnt)
    {
        status_t status = SLM_SUCCESS;

        if (isSem)
            sem_wait (&lock);

        if (filled < cap)
        {
            tail = plusOne(tail);
            arr[tail] = elmnt;
            filled++;
            status = SLM_SUCCESS;
        }
        else
        {
            status = ISFULL;
        }

        if (isSem)
            sem_post (&lock);

        return status;
    }

    status_t pop()
    {
        status_t status;

        if (isSem)
            sem_wait (&lock);

        if (filled > 0)
        {
            head = plusOne(head);
            filled--;
            status = SLM_SUCCESS;
        }
        else
        {
            status = ISEMPTY;
        }

        if (isSem)
            sem_post (&lock);

        return status;
    }

    BOOL isEmpty()
    {
        BOOL res = false;

        if (isSem)
            sem_wait (&lock);

        res = (filled == 0);

        if (isSem)
            sem_post (&lock);

        return res;
    }

    BOOL isFull()
    {
        BOOL res = false;

        if (isSem)
            sem_wait (&lock);

        res = (filled == cap);

        if (isSem)
            sem_post (&lock);

        return res;
    }

    T front()
    {
        T res;

        if (isSem)
            sem_wait (&lock);

        if (filled > 0)
            res = arr[head];
        else
            res = 0;

        if (isSem)
            sem_post (&lock);

        return res;
    }

    T end()
    {
        T res;

        if (isSem)
            sem_wait (&lock);

        if (filled > 0)
            res = arr[tail];
        else
            res = 0;

        if (isSem)
            sem_post (&lock);

        return res;
    }

    int_t size()
    {
        int_t res = 0;

        if (isSem)
            sem_wait (&lock);

        res = filled;

        if (isSem)
            sem_post (&lock);

        return res;
    }


};