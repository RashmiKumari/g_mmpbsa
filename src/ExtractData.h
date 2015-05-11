/*
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
 *
 */


#ifndef EXTRACTDATA_H_
#define EXTRACTDATA_H_

typedef int bool;

#ifndef FALSE
#  define FALSE   0
#endif
#ifndef TRUE
#  define TRUE    1
#endif

#ifndef STRLEN
#define STRLEN 256
#endif

double* extract_coulmn_double(char *str, int col_min, int col_max);
/*
 return pointer of array of size [(col_max-col_min)+1]
 array contains floating point numbers range from from col_min to col_max of the given line/string.
 For example;
 ---------------------------------------
 	double *out;
 	char *line = "    24 U-A       0.07      0.02      0.49     -2.48     -7.63     -0.38";
 	char *input=NULL;
 	input = strdup(line);
	out = extract_coulmn_double(input, 3, 5);
	printf("%15.3f%15.3f%15.3f\n",out[0], out[1], out[2]);
-----------------------------------------
	OUTPUT $          0.070          0.020          0.490
 */

int* extract_coulmn_integer(char *str, int col_min, int col_max);
/*
 Similar to extract_coulmn_double but extract integer data.
 */

char** get_all_lines(FILE *fp, int *num_line);
/*
 Return lines from the input file stream.
 Each line/string end with '\0' null terminator.
 num_line is output and store the total number of line in the given file.
 For example:
 -----------------------------------------
 	FILE *fp=NULL;
	char **lines=NULL;
	int num_line = 0, i =0;
	fp = fopen("any_file.dat", "r");
	lines = get_all_lines(fp, &num_line);			//Pointer of num_line is passed

	printf("Total number of lines in file: %d",num_line);
	for (i=0;i<num_line;i++){
		printf("%s",lines[i]);
	}
-------------------------------------------
 */

char* get_line(FILE *fp);
/*
 Return one line from the file
 This function is used in get_all_lines.
 Each line/string end with '\0' null terminator.
 */

void remove_leading_white_space(char *str);
/*
 Remove any leading white space from the input line/string
 and replace input string with new line/string.
 */

char** split_by_space(char *str,  int *num);
/*
 Split input line with the white spaces
 and return the list of words in the given line


 "num" is total number of string in the array. It
 could be NULL also.

 */

char** split_by_char(char *str, char *delimeter, int *num);
/*
 Similar to split sub-routine in perl
 Split input line with the given character
 and return the list of words in the given line


 "num" is total number of string in the array. It
 could be NULL also.


 For example:
 -------------------------------
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
 ---------------------------------
 OUTPUT $..51_ ## ..63_
 OUTPUT $String Type: 51
 OUTPUT $Integer Type: 51
 */

char* extract_digits(char *str);
/*
 Extract integer embedded in the given word string.
 After extracting the digits, one can directly convert it to integer type bye atoi() function
 For example:
 -------------------------------
 char **split_data=NULL;
 char *line = "10   (0.034) ...1>-:..51_:[.RC]C-----G[.RG]:..63_:-<...1 (0.040)     |"
 char *input = NULL:
 char *only_digit=NULL;

 input = strdup(line);
 split_data = split_by_char(input, ":");				//Splitting with semi-colon.
 only_digits = extract_digits(split_data[1]);			//Extracting digits from "..51_"
 printf("String Type: %s\n",only_digits);
 printf("Integer Type: %d\n",atoi(only_digits));
 ---------------------------------
 OUTPUT $ 51

 */

bool is_first_numerics(char *str);

char** get_block_lines(FILE *fp, char *delimeter, int *num_line);

#endif /* STRING_MANIP_H_ */
