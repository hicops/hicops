/*
 * Copyright (C) 2021  Muhammad Haseeb and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is licensed under the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * See full license terms here: http://creativecommons.org/licenses/by-nc-nd/4.0/
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
