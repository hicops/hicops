/*
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

#pragma once

/* Error Types and Codes Here */

#define SLM_SUCCESS               0

#define SLM_ERROR_BASE           -3000

/* Proper errors that must be thrown when something wrong happens
 * Define more errors as the development goes
 */
#define ERR_INVLD_MEMORY         (SLM_ERROR_BASE -  1)
#define ERR_INVLD_PARAM          (SLM_ERROR_BASE -  2)
#define ERR_INVLD_INDEX          (SLM_ERROR_BASE -  3)
#define ERR_INVLD_SIZE           (SLM_ERROR_BASE -  4)
#define ERR_SLM_ALREADY_INIT     (SLM_ERROR_BASE -  5)
#define ERR_FILE_NOT_FOUND       (SLM_ERROR_BASE -  6)
#define ERR_BAD_MEM_ALLOC        (SLM_ERROR_BASE -  7)
#define ERR_INVLD_PTR            (SLM_ERROR_BASE -  8)
#define ERR_INVLD_SORT           (SLM_ERROR_BASE -  9)
#define ERR_INVLD_MOD            (SLM_ERROR_BASE - 10)
#define ERR_INVLD_NODE_RANK      (SLM_ERROR_BASE - 11)
#define ERR_NOT_ENOUGH_DATA      (SLM_ERROR_BASE - 12)
#define ERR_FILE_ERROR           (SLM_ERROR_BASE - 13)


#define ENDSIGNAL                (SLM_ERROR_BASE - 100)
