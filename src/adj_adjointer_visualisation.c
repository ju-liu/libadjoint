#include "libadjoint/adj_adjointer_visualisation.h"

void adj_html_css(FILE* fp) {
	fprintf(fp, "<head>\n"
				"<style type=\"text/css\">\n"
			"table.equations\n"
			"{ font-family: monospace;\n"
			"  font-weight: normal;\n"
			"  font-size: 11px;\n"
			"  color: #404040;\n"
			"  table-layout:fixed;\n"
			"  background-color: #fafafa;\n"
			"  border: 1px #6699CC solid;\n"
			"  border-collapse: collapse;\n"
			"  border-spacing: 0px;\n"
			"  margin-top: 0px;}\n"

			"table.equations td\n"
			"{  border-bottom: 1px #6699CC;\n"
			"	font-weight: normal;\n"
			"	font-size: 11px;\n"
			"	color: #404040;\n"
			"	background-color: white;\n"
			"	text-align: left;\n"
			"	padding-left: 3px;}\n"

			"table.equations tr.new_iteration\n"
			"{ border-top: 2px solid #6699CC;}\n"

			"</style>\n"
			"</head>\n"
			);
}

void adj_html_header(FILE *fp) {
	adj_html_css(fp);
	fprintf(fp, "<html>\n"
				"<body>\n"
			);
}

void adj_html_footer(FILE *fp) {
	fprintf(fp, "</body>\n"
				"</html>\n"
			);
}

void adj_html_table_begin(FILE* fp){
	fprintf(fp, "<table border=\"1px\" class=\"equations\">\n");
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

int adj_html_find_column_index(adj_adjointer* adjointer, adj_variable* variable, int* col, int backward)
{
	int i;
	adj_equation* adj_eqn;
	for (i=0; i < adjointer->nequations; i++)
	{
		adj_eqn = &adjointer->equations[i];
		if (adj_variable_equal(&adj_eqn->variable, variable, 1))
		{
			if (backward)
				*col = adjointer->nequations-i-1;
			else
				*col = i;
			return ADJ_ERR_OK;
		}
	}
	strncpy(adj_error_msg, "Variable not found.", ADJ_ERROR_MSG_BUF);
	return ADJ_ERR_INVALID_INPUTS;
}

/* Writes a html row containing the adjoint variables into fp */
void adj_html_vars(FILE* fp, adj_adjointer* adjointer, int type, int backward)
{
	int i;
	char adj_name[ADJ_NAME_LEN];
	adj_variable adj_var;
	fprintf(fp, "<tr>\n");
	for (i = 0; i < adjointer->nequations; i++)
	{
		if (backward)
			adj_var = adjointer->equations[adjointer->nequations-i-1].variable;
		else
			adj_var = adjointer->equations[i].variable;
		adj_var.type = type;
		adj_variable_str(adj_var, adj_name, ADJ_NAME_LEN);
		fprintf(fp, "<th width=\"500px\">%s:%d:%d</th>\n", adj_name, adj_var.timestep, adj_var.iteration);
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
		ierr = adj_html_find_column_index(adjointer, &adj_eqn.targets[i], &col, ADJ_FALSE);
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
	/* Write it to file */
	adj_html_write_row(fp, row, desc, adjointer->nequations);


	/* Tidy up */
	for (i = 0; i < adjointer->nequations; ++i)
	{
		free(row[i]);
		free(desc[i]);
	}
	return ADJ_ERR_OK;
}

/* Writes a html row containing the supplied adjoint equation into fp */
int adj_html_adjoint_eqn(FILE* fp, adj_adjointer* adjointer, adj_equation fwd_eqn){
	int i, j;
	char* row[adjointer->nequations];
	char* desc[adjointer->nequations];
	char buf[ADJ_NAME_LEN];
	int col, ierr;
	adj_variable fwd_var;
	adj_variable_data *fwd_data;


	/* Allocate the strings for this row */
	for (i = 0; i < adjointer->nequations; ++i)
	{
		row[i] = malloc(ADJ_NAME_LEN*sizeof(char));
		row[i][0]='\0';
		desc[i] = malloc(ADJ_NAME_LEN*sizeof(char));
		desc[i][0]='\0';
	}

	fwd_var = fwd_eqn.variable;
	ierr = adj_find_variable_data(&(adjointer->varhash), &fwd_var, &fwd_data);
	assert(ierr == ADJ_ERR_OK);

	/* Each targeting equation corresponds to one column in the considered row */
	for (i = 0; i < fwd_data->ntargeting_equations; i++)
	{
	    adj_equation other_fwd_eqn;

	    other_fwd_eqn = adjointer->equations[fwd_data->targeting_equations[i]];

	    /* Find the index in other_fwd_eqn of the block fwd_var is multiplied with */
	    for (j=0; j<other_fwd_eqn.nblocks; j++)
	    {
	    	if (adj_variable_equal(&other_fwd_eqn.targets[j], &fwd_var, 1))
	    		break;
	    }
	    assert(j!=other_fwd_eqn.nblocks);

	    /* find the column in which this blocks belongs */
	    ierr = adj_html_find_column_index(adjointer, &other_fwd_eqn.variable, &col, ADJ_TRUE);
	    if (ierr != ADJ_ERR_OK)
	    	return ierr;

	    //other_adj_var.type = ADJ_ADJOINT;


		/* Fill in the data */
		strncpy(row[col], other_fwd_eqn.blocks[j].name, ADJ_NAME_LEN);

		strncpy(desc[col], other_fwd_eqn.blocks[j].name, ADJ_NAME_LEN);
		strncat(desc[col], "\n\nTargets: ", ADJ_NAME_LEN);
		strncat(desc[col], other_fwd_eqn.targets[j].name, ADJ_NAME_LEN);
		strncat(desc[col], ":", ADJ_NAME_LEN);
		snprintf(buf, ADJ_NAME_LEN, "%d", other_fwd_eqn.targets[j].timestep);
		strncat(desc[col], buf, ADJ_NAME_LEN);
		strncat(desc[col], ":", ADJ_NAME_LEN);
		snprintf(buf, ADJ_NAME_LEN, "%d", other_fwd_eqn.targets[j].iteration);
		strncat(desc[col], buf, ADJ_NAME_LEN);
		strncat(desc[col], "\nCoefficient: ", ADJ_NAME_LEN);
		snprintf(buf, ADJ_NAME_LEN, "%f", other_fwd_eqn.blocks[j].coefficient);
		strncat(desc[col], buf, ADJ_NAME_LEN);
		if (other_fwd_eqn.blocks[j].has_nonlinear_block)
		{
			strncat(desc[col], "\nNonlinear Block: ", ADJ_NAME_LEN);
			strncat(desc[col], other_fwd_eqn.blocks[j].nonlinear_block.name, ADJ_NAME_LEN);
		}

	}

	/* Write it to file */
	adj_html_write_row(fp, row, desc, adjointer->nequations);

	/* Tidy up */
	for (i = 0; i < adjointer->nequations; ++i)
	{
		free(row[i]);
		free(desc[i]);
	}
	return ADJ_ERR_OK;
}

int adj_html_adjoint_system(adj_adjointer* adjointer, char* filename)
{
	FILE *fp;
	int i, ierr, timestep, iteration;
	adj_equation adj_eqn;

	fp=fopen(filename, "w+");
	if (fp==NULL)
	{
		strncpy(adj_error_msg, "File could not be opened.", ADJ_ERROR_MSG_BUF);
		return ADJ_ERR_INVALID_INPUTS;
	}

	adj_html_header(fp);

	fprintf(fp, "<h1>Registered equations</h1>\n");
	fprintf(fp, "Adjoint system\n");
	if (adjointer->nequations==0)
		return ADJ_ERR_OK;
	timestep = adjointer->equations[0].variable.timestep-1;
	iteration = adjointer->equations[0].variable.iteration - 1 ;

	/* The adjoint equation runs backward in time */
	for (i=adjointer->nequations-1; i>=0; i--)
	{
		adj_eqn = adjointer->equations[i];
		if (timestep != adj_eqn.variable.timestep) {
			/* New timestep, create a new table */
			if (i!=adjointer->nequations-1)
				adj_html_table_end(fp);
			fprintf(fp, "<h2>Timestep %d</h2>\n", adj_eqn.variable.timestep);
			adj_html_table_begin(fp);
			adj_html_vars(fp, adjointer, ADJ_ADJOINT, ADJ_TRUE);
			fprintf(fp, "<tr class=\"new_iteration\">\n");
		}
		else {
			if (iteration != adj_eqn.variable.iteration)
				fprintf(fp, "<tr class=\"new_iteration\">\n");
			else
				fprintf(fp, "<tr>\n");
		}

		ierr = adj_html_adjoint_eqn(fp, adjointer, adj_eqn);
		if (ierr != ADJ_ERR_OK)
		{
			fclose(fp);
			return ierr;
		}
		fprintf(fp, "</tr>\n");
		timestep = adj_eqn.variable.timestep;
		iteration = adj_eqn.variable.iteration;
	}
	adj_html_table_end(fp);

	adj_html_footer(fp);
	fclose(fp);
	return ADJ_ERR_OK;
}


int adj_html_forward_system(adj_adjointer* adjointer, char* filename)
{
	FILE *fp;
	int i, ierr, timestep, iteration;
	adj_equation adj_eqn;

	fp=fopen(filename, "w+");
	if (fp==NULL)
	{
		strncpy(adj_error_msg, "File could not be opened.", ADJ_ERROR_MSG_BUF);
		return ADJ_ERR_INVALID_INPUTS;
	}

	adj_html_header(fp);

	fprintf(fp, "<h1>Registered equations</h1>\n");
	fprintf(fp, "Forward system\n");
	if (adjointer->nequations==0)
		return ADJ_ERR_OK;
	timestep = adjointer->equations[0].variable.timestep - 1;
	iteration = adjointer->equations[0].variable.iteration - 1 ;

	for (i=0; i<adjointer->nequations; i++)
	{
		adj_eqn = adjointer->equations[i];
		if (timestep != adj_eqn.variable.timestep) {
			/* New timestep, create a new table */
			if (i!=0)
				adj_html_table_end(fp);
			fprintf(fp, "<h2>Timestep %d</h2>\n", adj_eqn.variable.timestep);
			adj_html_table_begin(fp);
			adj_html_vars(fp, adjointer, ADJ_FORWARD, ADJ_FALSE);
			fprintf(fp, "<tr class=\"new_iteration\">\n");
		}
		else {
			if (iteration != adj_eqn.variable.iteration)
				fprintf(fp, "<tr class=\"new_iteration\">\n");
			else
				fprintf(fp, "<tr>\n");
		}

		ierr = adj_html_eqn(fp, adjointer, adj_eqn);
		if (ierr != ADJ_ERR_OK)
		{
			fclose(fp);
			return ierr;
		}
		fprintf(fp, "</tr>\n");
		timestep = adj_eqn.variable.timestep;
		iteration = adj_eqn.variable.iteration;
	}
	adj_html_table_end(fp);

	adj_html_footer(fp);
	fclose(fp);
	return ADJ_ERR_OK;
}

int adj_adjointer_to_html(adj_adjointer* adjointer, char* filename, int type)
{
	if (type == ADJ_FORWARD)
		return adj_html_forward_system(adjointer, filename);
	else if(type == ADJ_ADJOINT)
		return adj_html_adjoint_system(adjointer, filename);
	else
		return ADJ_ERR_INVALID_INPUTS;
}
