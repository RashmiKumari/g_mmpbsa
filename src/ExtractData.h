/*!
 * @defgroup Ext_apbs  Routins to extract data from APBS output file
 */

/*!
 *
 *  @file    ExtractData.h
 *  @author  Rajendra Kumar
 *  @brief   Header file for routins to extract data from APBS output files
 *  @ingroup Ext_apbs
 *  @attention
 *  @verbatim
 *
 * This file is part of g_mmpbsa.
 *
 * Authors: Rajendra Kumar
 *
 * Copyright (C) 2013, 2014, 2015 Rajendra Kumar
 *
 * g_mmpbsa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * g_mmpbsa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with g_mmpbsa.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */


#ifndef EXTRACTDATA_H_
#define EXTRACTDATA_H_

/**
 * @brief Define `bool` variable with `int`
 * @ingroup Ext_apbs
 */
typedef int bool;

/**
 * @brief Define bool variable `FALSE`
 * @ingroup Ext_apbs
 */
#ifndef FALSE
#  define FALSE   0
#endif

/**
 * @brief Define bool variable `TRUE`
 * @ingroup Ext_apbs
 */
#ifndef TRUE
#  define TRUE    1
#endif

/**
 * @brief Define maximum string length of 256
 * @ingroup Ext_apbs
 */
#ifndef STRLEN
#define STRLEN 256
#endif

/**
 * @brief To get array of `double` from a string, which contains values in column
 * @ingroup Ext_apbs
 * @param[in] str Input string
 * @param[in] col_min column number from which values will be extracted
 * @param[in] col_max column number up to which values will be extracted
 * @returns pointer to a array of double of size `[(col_max-col_min)+1]`. Array contains floating point numbers range from from `col_min` to `col_max` of the given line/string.

 <b> For example: </b>
\code{.c}
    double *out;
    char *line = "    24 U-A       0.07      0.02      0.49     -2.48     -7.63     -0.38";
    char *input=NULL;
    input = strdup(line);
    out = extract_coulmn_double(input, 3, 5);
    printf("%15.3f%15.3f%15.3f\n",out[0], out[1], out[2]);
\endcode

\b OUTPUT: `0.070 0.020 0.490`

 */
double* extract_coulmn_double(char *str, int col_min, int col_max );

/**
 * @ingroup Ext_apbs
 * @brief Similar to extract_coulmn_double() but extract integer data.
 * @see extract_coulmn_double
 */
int* extract_coulmn_integer(char *str, int col_min, int col_max);

/**
 * @ingroup Ext_apbs
 * @brief Return array of strings/lines from the input file stream. Each line/string end with '\0' null terminator.
 * @param[in] fp input file stream
 * @param[out] num_line Output number of lines/strings read from input file stream (pass by reference). It is output and store the total number of line in the given file.
 * @returns arrays of strings of size `num_line`

 <b>For example:</b>
\code{.c}
 	FILE *fp=NULL;
	char **lines=NULL;
	int num_line = 0, i =0;
	fp = fopen("any_file.dat", "r");
	lines = get_all_lines(fp, &num_line);			//Pointer of num_line is passed

	printf("Total number of lines in file: %d",num_line);
	for (i=0;i<num_line;i++){
		printf("%s",lines[i]);
	}
\endcode
*/
char** get_all_lines(FILE *fp, int *num_line);

/**
 * @ingroup Ext_apbs
 * @brief To get a line or string from the input file stream. This function is used in get_all_lines()
 * @param[in] fp input file stream
 * @returns line/string end with <code>'\0'</code> null terminator.
 */
char* get_line(FILE *fp);

/**
 * @ingroup Ext_apbs
 * @brief Remove any leading white space from the input line/string and replace input string with new line/string.
 * @param[in] str Input string/line
 * @returns None
 */
void remove_leading_white_space(char *str);

/**
 * @ingroup Ext_apbs
 * @brief To split input line/string with the white space delimiter
 * @param[in] str input string/line
 * @param[out] num number (output) of words in string or length of output array. It could be `NULL` also.
 * @returns the array of words/strings of length `num`
 */
char** split_by_space(char *str,  int *num);

/**
 * @ingroup Ext_apbs
 * @brief Similar to split sub-routine in perl. Split input line with the given input character and return the array of words.
 * @param[in] str Input string/line.
 * @param[in] delimeter character at which line will be spllited.
 * @param[out] num Total number of string/words in the returned array. It could be `NULL` also.
 * @returns A array of words or strings of length `num`.

 <b>For example:</b>
\code{.c}
 char **split_data=NULL;
 char *line = "10   (0.034) ...1>-:..51_:[.RC]C-----G[.RG]:..63_:-<...1 (0.040)     |";
 char *input=NULL;
 input = strdup(line);
 char *only_digits=NULL;

 split_data = split_by_char(input, ":");				//Splitting with semi-colon.
 printf("%s ## %s\n",split_data[1], split_data[3]);
 only_digits = extract_digits(split_data[1]);			//Extracting digits from "..51_"
 printf("String Type: %s\n",only_digits);
 printf("Integer Type: %d\n",atoi(only_digits));
\endcode

\b OUTPUT:
\code
..51_ ## ..63_
String Type: 51
Integer Type: 51
\endcode
 */
char** split_by_char(char *str, char *delimeter, int *num);

/**
 * @brief Extract integer present within alphanumeric word/string.
 * @ingroup Ext_apbs
 * @param[in] str Input word/string
 * @returns a integer of string type. one can directly convert it to integer type by `atoi()`

 <b> For example:</b>
 \code
 char **split_data=NULL;
 char *line = "10   (0.034) ...1>-:..51_:[.RC]C-----G[.RG]:..63_:-<...1 (0.040)     |"
 char *input = NULL:
 char *only_digit=NULL;

 input = strdup(line);
 split_data = split_by_char(input, ":");				//Splitting with semi-colon.
 only_digits = extract_digits(split_data[1]);			//Extracting digits from "..51_"
 printf("String Type: %s\n",only_digits);
 printf("Integer Type: %d\n",atoi(only_digits));
 \endcode

 \b OUTPUT: 51
 */
char* extract_digits(char *str);

/**
* @brief To check if first charecter of line is digit
* @ingroup Ext_apbs
* @param[in] str Input string or line
* @returns `TRUE` or `FALSE`
*/
bool is_first_numerics(char *str);

/**
 * @ingroup Ext_apbs
 * @brief To get strings/lines from the input file stream till `delimiter` is found. A specific section or block can be read from file.
 * @param[in] fp input file stream
 * @param[in] delimeter A word or string till file should be read
 * @param[out] num_line Output number of lines/strings read from input file stream (pass by reference). It is output and store the total number of line in the given section/block.
 * @returns arrays of strings of size `num_line`

  This function is very useful for file which contains a identifier for each section or block. For example, a `"Time = "` identifier shows that the next section is output from the next frame.
  This function would return all the lines corresponding to the frame.
*/
char** get_block_lines(FILE *fp, char *delimeter, int *num_line);

#endif
