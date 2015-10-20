/*
 * This file is part of g_mmpbsa.
 *
 * Authors: Rajendra Kumar
 *
 * Copyright (C) 2013-2015 Rajendra Kumar
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "ExtractData.h"

int* extract_coulmn_integer(char *str, int col_min, int col_max)	{
	/*
	 Similar to extract_coulmn_double but extract integer data.
	 */

	int *data=NULL;
	char *buffer=NULL, **str_data=NULL;
	int size = (col_max-col_min)+1;
	int i=0,n=0;

	buffer = strdup(str);
	remove_leading_white_space(buffer);
	str_data = split_by_space(buffer, NULL);

	data = (int *) malloc (sizeof(int)*size);
	for(i=(col_min-1);i<col_max;i++){
		data[n] = atoi(str_data[i]);
		n++;
	}
	return data;
}


double* extract_coulmn_double(char *str, int col_min, int col_max)	{
	/*
	 return pointer to array of size [(col_max-col_min)+1]
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

	double *data=NULL;
	char *buffer=NULL, **str_data=NULL;
	int size = (col_max-col_min)+1;
	int i=0,n=0;

	buffer = strdup(str);
	remove_leading_white_space(buffer);
	str_data = split_by_space(buffer, NULL);

	data = (double *) malloc (sizeof(double)*size);
	for(i=(col_min-1);i<col_max;i++){
		data[n] = strtof(str_data[i],NULL);
		n++;
	}
	free(buffer);
	free(str_data);
	return data;
}

bool is_first_numerics(char *str)	{
	char *buffer;
	bool inumber=FALSE;
	buffer = strdup(str);
	remove_leading_white_space(buffer);
	if(isdigit(buffer[0]))
		inumber=TRUE;
	if((buffer[0]=='-') && (isdigit(buffer[1])))
			inumber=TRUE;
	free(buffer);
	return inumber;
}

char* extract_digits(char *str)	{
	/*
	 Extract integer embedded in the given word string.
	 After extracting the digits, one can directly convert it to integer type by atoi() function
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
	 OUTPUT $   51

	 */

	char *final=NULL;
	int i = 0, n=0;
	int len = strlen(str)+1;

	final = (char *) malloc (sizeof(char)*len);

	do {

		if isdigit(str[i])	{
			final[n] = str[i];
			n++;
		}

		i++;
	} while (i<len);

	final[n] = '\0';

	return final;
}

char** split_by_char(char *str, char *delimeter, int *num)	{
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

	char **final = NULL;
	char *token;
	char *input = str;
	int n = 0;

	input = strdup(str);
	final = (char **) malloc (sizeof(char*));

	if(input!=NULL)	{

		token = strtok(input, delimeter);

		while(token != NULL)	{
			final  = (char**) realloc (final, (sizeof(char*)*(n+1)));

			final[n] = token;
			n++;
			token = strtok(NULL, delimeter);
		}
	}

	if (num!=NULL)
		*num = n;

	return final;
}


char** split_by_space(char *str, int *num)	{
	/*
	 Split input line with the white spaces
	 and return the list of words in the given line.

	 "num" is total number of string in the array. It
	 could be NULL also.
	 */

	char **final = NULL;
	char *token = NULL;
	char *input = NULL;
	int n = 0;

	input = strdup(str);
	final = (char **) malloc (sizeof(char*));

	if(input!=NULL)	{

		token = strtok(input, " \t\n\v\f\r");

		while(token != NULL)	{
			final  = (char**) realloc (final, (sizeof(char*)*(n+1)));

			final[n] = token;
			n++;
			token = strtok(NULL, " \t\n\v\f\r");
		}
	}

	if (num!=NULL)
		*num = n;

	return final;
}

void remove_leading_white_space(char *str)	{
	/*
	 Remove any leading white space from the input line/string
	 and replace input string with new line/string.
	 */

	char *final=NULL;
	int i = 0, nspace=0;
	int len = strlen(str)+1;

	final = (char *) malloc (sizeof(char)*len);

	for(i=0; i<=len; i++)	{

		if isspace(str[i])	{
			nspace += 1;
			continue;
		}
		else
			break;
	}

	for (i=nspace; i<=len; i++)	{

		final[i-nspace] = str[i];

		if (str[i]=='\0')
			break;
	}

	strcpy(str,final);

	free(final);
}

char* get_line(FILE *fp)	{
	/*
	 Return one line from the file
	 This function is used in get_all_lines.
	 Each line/string end with '\0' null terminator.
	 */

	char c, *str;
	char *final;
	int count = 0, mult=1, len;

	str = (char *) malloc (sizeof(char)*STRLEN);

	do {

		//Getting character from file
		c = fgetc(fp);

		if ((mult*STRLEN)<count)	{
			mult += 1;
			str  = (char*) realloc (str, (sizeof(char)*STRLEN*mult));
		}

		//Checking end of file
		if (c==EOF)	{
			str[count] = '\n';
			str[count+1] = '\0';
			count += 2;
			break;
		}

		//Checking end of the line
		if(c=='\n') 	{
			str[count] = c;
			str[count+1] = '\0';
			count += 2;
			break;
		}

		str[count] = c;

		count++;

	} while (c != EOF);

	len = strlen(str)+1;
	final = (char*) malloc (sizeof(char)*len);
	strcpy(final, str);
	free(str);

	return final;
}

char** get_all_lines(FILE *fp, int *num_line)	{
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

	char *buffer=NULL;
	int count = 0;
	char **data = 0;

	data = (char **) malloc (sizeof(char*));
	do {
			buffer = get_line(fp);

			if (buffer==NULL)
				continue;

			if(count==0)	{
				data[count] = buffer;

			}
			else	{
				data  = (char**) realloc (data, (sizeof(char*)*(count+1)));
				data[count] = buffer;
			}

			count++;

	} while(!feof(fp));


	*num_line = (count-1);
	return data;
}

char** get_block_lines(FILE *fp, char *delimeter, int *num_line)	{

	char *buffer=NULL;
	int count = 0;
	char **data = 0;

	data = (char **) malloc (sizeof(char*));
	while(1) {

		if(feof(fp))
			break;

			buffer = get_line(fp);

			if (buffer==NULL)
				continue;

			if(strstr(buffer,delimeter)!=NULL)
				break;

			if(count==0)	{
				data[count] = buffer;

			}
			else	{
				data  = (char**) realloc (data, (sizeof(char*)*(count+1)));
				data[count] = buffer;
			}

			count++;
	}

	if(buffer!=NULL){
		data  = (char**) realloc (data, (sizeof(char*)*(count+1)));
		data[count] = buffer;
		count++;
	}

	*num_line = (count-1);
	return data;
}
