/*
 * This file is part of PCDSFrame software
 * Copyright (C) 2019  Muhammad Haseeb, Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef LWQUEUE_H_
#define LWQUEUE_H_

#include "common.h"
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
    INT filled;
    INT    cap;
    INT   head;
    INT   tail;
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

    lwqueue(INT dcap)
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

    lwqueue(INT dcap, BOOL sem)
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

    inline INT plusOne(INT idx)
    {
        return ((idx + 1) % cap);
    }

    inline INT minusOne(INT idx)
    {
        return ((idx - 1) % cap);
    }

    STATUS push(T elmnt)
    {
        STATUS status = SLM_SUCCESS;

        if (isSem)
        {
            sem_wait (&lock);
        }

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
        {
            sem_post (&lock);
        }

        return status;
    }

    STATUS pop()
    {
        STATUS status;

        if (isSem)
        {
            sem_wait (&lock);
        }

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
        {
            sem_post (&lock);
        }

        return status;
    }

    BOOL isEmpty()
    {
        BOOL res = false;

        if (isSem)
        {
            sem_wait (&lock);
        }

        res = (filled == 0);

        if (isSem)
        {
            sem_post (&lock);
        }

        return res;
    }

    BOOL isFull()
    {
        BOOL res = false;

        if (isSem)
        {
            sem_wait (&lock);
        }

        res = (filled == cap);

        if (isSem)
        {
            sem_post (&lock);
        }

        return res;
    }

    T front()
    {
        T res;// = (T) NULL;

        if (isSem)
        {
            sem_wait (&lock);
        }

        if (filled > 0)
        {
            res = arr[head];
        }

        if (isSem)
        {
            sem_post (&lock);
        }

        return res;
    }

    T end()
    {
        T res;// = (T) NULL;

        if (isSem)
        {
            sem_wait (&lock);
        }

        if (filled > 0)
        {
            res = arr[tail];
        }

        if (isSem)
        {
            sem_post (&lock);
        }

        return res;
    }

    INT size()
    {
        INT res = 0;

        if (isSem)
        {
            sem_wait (&lock);
        }

        res = filled;

        if (isSem)
        {
            sem_post (&lock);
        }

        return res;
    }


};

#endif /* LWQUEUE_H_ */
