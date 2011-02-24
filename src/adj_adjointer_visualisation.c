#include "libadjoint/adj_adjointer_visualisation.h"

void adj_html_table_begin(FILE* fp){
	fprintf(fp, "<table border=\"1px\" style=\"border-collapse: collapse; border: solid;\">\n");
}

void adj_html_table_end(FILE* fp){
	fprintf(fp, "</table>\n");
}

void adj_html_write_row(FILE* fp, char** strings, char** desc, int nb_strings)
{
	int i;
	for (i = 0; i < nb_strings; i++)
		fprintf(fp, "<td><div title=\"%s\">%s</div></td>\n", desc[i], strings[i]);
}

int adj_html_find_column_index(adj_adjointer* adjointer, adj_variable* variable, int* col)
{
	int i;
	adj_equation* adj_eqn;
	for (i=0; i < adjointer->nequations; i++)
	{
		adj_eqn = &adjointer->equations[i];
		if (adj_variable_equal(&adj_eqn->variable, variable, 1))
		{
			*col = i;
			return ADJ_ERR_OK;
		}
	}
	strncpy(adj_error_msg, "Variable not found.", ADJ_ERROR_MSG_BUF);
	return ADJ_ERR_INVALID_INPUTS;
}

/* Writes a html row containing the variables in the adjointer into fp */
void adj_html_vars(FILE* fp, adj_adjointer* adjointer)
{
	int i;
	adj_equation* adj_eqn;
	fprintf(fp, "<tr>\n");
	for (i=0; i < adjointer->nequations; i++)
	{
		adj_eqn = &adjointer->equations[i];
		fprintf(fp, "<th>%s:%d:%d</th>\n", adj_eqn->variable.name, adj_eqn->variable.timestep, adj_eqn->variable.iteration);
	}
	fprintf(fp, "</tr>\n");
}

/* Writes a html row containing the supplied equation into fp */
int adj_html_eqn(FILE* fp, adj_adjointer* adjointer, adj_equation adj_eqn){
	int i;
	char* row[adjointer->nequations];
	char* desc[adjointer->nequations];
	char buf[ADJ_NAME_LEN];
	int col, ierr;

	/* Allocate the strings for this row */
	for (i = 0; i < adjointer->nequations; ++i)
	{
		row[i] = malloc(ADJ_NAME_LEN*sizeof(char));
		row[i][0]='\0';
		desc[i] = malloc(ADJ_NAME_LEN*sizeof(char));
		desc[i][0]='\0';
	}

	/* Fill in the data */
	for (i = 0; i < adj_eqn.nblocks; i++)
	{
		ierr = adj_html_find_column_index(adjointer, &adj_eqn.targets[i], &col);
		if (ierr != ADJ_ERR_OK)
			return ierr;

		strncpy(row[col], adj_eqn.blocks[i].name, ADJ_NAME_LEN);
		/* Fill in the description */
		strncpy(desc[col], adj_eqn.blocks[i].name, ADJ_NAME_LEN);

		strncat(desc[col], "\n\nTargets: ", ADJ_NAME_LEN);
		strncat(desc[col], adj_eqn.targets[i].name, ADJ_NAME_LEN);
		strncat(desc[col], ":", ADJ_NAME_LEN);
		snprintf(buf, ADJ_NAME_LEN, "%d", adj_eqn.targets[i].timestep);
		strncat(desc[col], buf, ADJ_NAME_LEN);
		strncat(desc[col], ":", ADJ_NAME_LEN);
		snprintf(buf, ADJ_NAME_LEN, "%d", adj_eqn.targets[i].iteration);
		strncat(desc[col], buf, ADJ_NAME_LEN);

		strncat(desc[col], "\nCoefficient: ", ADJ_NAME_LEN);
		snprintf(buf, ADJ_NAME_LEN, "%f", adj_eqn.blocks[i].coefficient);
		strncat(desc[col], buf, ADJ_NAME_LEN);

		if (adj_eqn.blocks[i].has_nonlinear_block)
		{
			strncat(desc[col], "\nNonlinear Block: ", ADJ_NAME_LEN);
			strncat(desc[col], adj_eqn.blocks[i].nonlinear_block.name, ADJ_NAME_LEN);
		}
	}
	/* And write it to file */
	adj_html_write_row(fp, row, desc, adjointer->nequations);


	/* Tidy up */
	for (i = 0; i < adjointer->nequations; ++i)
	{
		free(row[i]);
		free(desc[i]);
	}
	return ADJ_ERR_OK;

}

int adj_adjointer_to_html(adj_adjointer* adjointer, char* filename)
{
	FILE *fp;
	int i, ierr;
	adj_equation adj_eqn;

	fp=fopen(filename, "w+");
	if (fp==NULL)
	{
		strncpy(adj_error_msg, "File could not be opened.", ADJ_ERROR_MSG_BUF);
		return ADJ_ERR_INVALID_INPUTS;
	}

	fprintf(fp, "<h1>Registered equations</h1>");
	adj_html_table_begin(fp);
	adj_html_vars(fp, adjointer);
	for (i=0; i<adjointer->nequations; i++)
	{
		adj_eqn = adjointer->equations[i];
		fprintf(fp, "<tr>\n");
		ierr = adj_html_eqn(fp, adjointer, adj_eqn);
		if (ierr != ADJ_ERR_OK)
		{
			fclose(fp);
			return ierr;
		}
		fprintf(fp, "</tr>\n");
	}
	adj_html_table_end(fp);

	fclose(fp);
	return ADJ_ERR_OK;
}
