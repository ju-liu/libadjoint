#include "libadjoint/adj_adjointer_visualisation.h"

void adj_html_css(FILE* fp)
{
  fprintf(fp, "<head>\n"
        "<style type=\"text/css\">\n"
      "table.equations\n"
      "{ font-family: monospace;\n"
      "  font-weight: normal;\n"
      "  font-size: 11px;\n"
      "  color: #404040;\n"
      "  table-layout:fixed;\n"
      "  background-color: #fafafa;\n"
      "  border: 0px;\n"
      "  empty-cells: hide;\n"
      "  border-spacing: 0px;\n"
      "  margin-top: 0px;}\n"

      "table.equations td\n"
      "{ font-weight: normal;\n"
      "  font-size: 11px;\n"
      "  color: #404040;\n"
      "  background-color: lightgray;\n"
      "  text-align: left;\n"
      "  padding-left: 3px;\n"
	  "  border:1px solid black;\n"
      "}\n"

      "table.equations td.new_timestep\n"
      "{ border-top: 5px solid black;}\n"

      "table.equations td.new_iteration\n"
      "{ border-top: 2px solid black}\n"

      "table.equations td.diagonal\n"
      "{ background-color: lightgreen;}\n"

	  ".box_rotate {\n"
	  "	       -moz-transform: rotate(310deg);  /* FF3.5+ */\n"
	  "	         -o-transform: rotate(310deg);  /* Opera 10.5 */\n"
	  "	    -webkit-transform: rotate(310deg);  /* Saf3.1+, Chrome */\n"
	  "	               filter:  progid:DXImageTransform.Microsoft.BasicImage(rotation=-0.0698131701);  /* IE6,IE7 */\n"
	  " 	           -ms-filter: 'progid:DXImageTransform.Microsoft.BasicImage(rotation=-0.0698131701)'; /* IE8 */\n"
      "}"

	  ".tablecell {\n"
	  "     width:40px;\n "
	  "     height:40px;}\n"

	  ".headercell {"
	  "     width:40px;\n"
	  "     height:40px;\n"
	  "     float:left;\n"
	  "     padding-left:5px;\n"
	  "}"

      "</style>\n"
      "</head>\n"
      );
}

void adj_html_header(FILE *fp)
{
  adj_html_css(fp);
  fprintf(fp, "<html>\n"
        "<body>\n"
      );
}

void adj_html_footer(FILE *fp)
{
  fprintf(fp, "</body>\n"
        "</html>\n"
      );
}

void adj_html_table_begin(FILE* fp, char* args)
{
  fprintf(fp, "<table class=\"equations\" %s>\n", args);
}

void adj_html_table_end(FILE* fp)
{
  fprintf(fp, "</table>\n");
}

void adj_html_write_row(FILE* fp, char** strings, char** desc, int nb_strings, int diag_index, char* class)
{
  int i;
  for (i = 0; i < nb_strings; i++)
  {
    if (strlen(desc[i]))
      if (diag_index == i)
        fprintf(fp, "<td class=\"diagonal tablecell %s\"><div title=\"%s\">%s</div></td>\n", class, desc[i], strings[i]);
      else
        fprintf(fp, "<td class=\"tablecell %s\"><div title=\"%s\">%s</div></td>\n", class, desc[i], strings[i]);
    else
      fprintf(fp, "<td class=\"tablecell %s\"></td>\n", class);
  }
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
      return ADJ_OK;
    }
  }
  strncpy(adj_error_msg, "Variable not found.", ADJ_ERROR_MSG_BUF);
  return ADJ_ERR_INVALID_INPUTS;
}

void adj_html_print_statistics(FILE* fp, adj_adjointer* adjointer)
{
	fprintf(fp, "<h1>General information</h1>");
	fprintf(fp, "Number of timesteps: %i<br>", adjointer->ntimesteps);
	fprintf(fp, "Number of registered equations: %i<br>", adjointer->nequations);
}

void adj_html_print_callback_information(FILE* fp, adj_adjointer* adjointer)
{
	adj_op_callback* cb_ptr;
	adj_func_callback* func_cb_ptr;
	adj_func_deriv_callback* func_deriv_cb_ptr;
	fprintf(fp, "<h1>Callback information</h1>");

	fprintf(fp, "<h2>Data callbacks</h2>");

	if (adjointer->callbacks.vec_duplicate == NULL)
		fprintf(fp, "<span style='color:red' title='Not registered'>vec_duplicate</span>");
	else
		fprintf(fp, "<span style='color:green'>vec_duplicate</span>");

	if (adjointer->callbacks.vec_axpy == NULL)
		fprintf(fp, ", <span style='color:red' title='Not registered'>vec_axpy</span>");
	else
		fprintf(fp, ", <span style='color:green'>vec_axpy</span>");

	if (adjointer->callbacks.vec_destroy == NULL)
		fprintf(fp, ", <span style='color:red' title='Not registered'>vec_destroy</span>");
	else
		fprintf(fp, ", <span style='color:green'>vec_destroy</span>");

	if (adjointer->callbacks.vec_set_values == NULL)
		fprintf(fp, ", <span style='color:red' title='Not registered'>vec_set_values</span>");
	else
		fprintf(fp, ", <span style='color:green'>vec_set_values</span>");

	if (adjointer->callbacks.vec_get_size == NULL)
		fprintf(fp, ", <span style='color:red' title='Not registered'>vec_get_size</span>");
	else
		fprintf(fp, ", <span style='color:green'>vec_get_size</span>");

	if (adjointer->callbacks.vec_divide == NULL)
		fprintf(fp, ", <span style='color:red' title='Not registered'>vec_divide</span>");
	else
		fprintf(fp, ", <span style='color:green'>vec_divide</span>");

	if (adjointer->callbacks.vec_get_norm == NULL)
		fprintf(fp, ", <span style='color:red' title='Not registered'>vec_get_norm</span>");
	else
		fprintf(fp, ", <span style='color:green'>vec_get_norm</span>");

	if (adjointer->callbacks.vec_dot_product == NULL)
		fprintf(fp, ", <span style='color:red' title='Not registered'>vec_dot_product</span>");
	else
		fprintf(fp, ", <span style='color:green'>vec_dot_product</span>");

	if (adjointer->callbacks.vec_set_random == NULL)
		fprintf(fp, ", <span style='color:red' title='Not registered'>vec_set_random</span>");
	else
		fprintf(fp, ", <span style='color:green'>vec_set_random</span>");

	if (adjointer->callbacks.mat_duplicate == NULL)
		fprintf(fp, ", <span style='color:red' title='Not registered'>mat_duplicate</span>");
	else
		fprintf(fp, ", <span style='color:green'>mat_duplicate</span>");

	if (adjointer->callbacks.mat_duplicate == NULL)
		fprintf(fp, ", <span style='color:red' title='Not registered'>mat_duplicate</span>");
	else
		fprintf(fp, ", <span style='color:green'>mat_duplicate</span>");

	if (adjointer->callbacks.mat_axpy == NULL)
		fprintf(fp, ", <span style='color:red' title='Not registered'>mat_axpy</span>");
	else
		fprintf(fp, ", <span style='color:green'>mat_axpy</span>");

	if (adjointer->callbacks.mat_destroy == NULL)
		fprintf(fp, ", <span style='color:red' title='Not registered'>mat_destroy</span>");
	else
		fprintf(fp, ", <span style='color:green'>mat_destroy</span>");


	fprintf(fp, "<h2>Block action callbacks</h2>");
	cb_ptr= adjointer->block_action_list.firstnode;
	while (cb_ptr != NULL) {
		fprintf(fp, "%s", cb_ptr->name);
		cb_ptr = cb_ptr->next;
	}

	fprintf(fp, "<h2>Block assembly callbacks</h2>");
	cb_ptr= adjointer->block_assembly_list.firstnode;
	while (cb_ptr != NULL) {
		fprintf(fp, "%s", cb_ptr->name);
		cb_ptr = cb_ptr->next;
	}

	fprintf(fp, "<h2>Nonlinear block action callbacks</h2>");
	cb_ptr= adjointer->nonlinear_action_list.firstnode;
	while (cb_ptr != NULL) {
		fprintf(fp, "%s", cb_ptr->name);
		cb_ptr = cb_ptr->next;
	}

	fprintf(fp, "<h2>Nonlinear derivative block action assembly callbacks</h2>");
	cb_ptr= adjointer->nonlinear_derivative_action_list.firstnode;
	while (cb_ptr != NULL) {
		fprintf(fp, "%s", cb_ptr->name);
		cb_ptr = cb_ptr->next;
	}

	fprintf(fp, "<h2>Nonlinear derivative block assembly assembly callbacks</h2>");
	cb_ptr= adjointer->nonlinear_derivative_assembly_list.firstnode;
	while (cb_ptr != NULL) {
		fprintf(fp, "%s", cb_ptr->name);
		cb_ptr = cb_ptr->next;
	}

	fprintf(fp, "<h2>Functional callbacks</h2>");
	func_cb_ptr= adjointer->functional_list.firstnode;
	while (func_cb_ptr != NULL) {
		fprintf(fp, "<span style='color:green'>%s</span>", func_cb_ptr->name);
		func_cb_ptr = func_cb_ptr->next;
	}

	fprintf(fp, "<h2>Functional derivative callbacks</h2>");
	func_deriv_cb_ptr= adjointer->functional_derivative_list.firstnode;
	while (func_deriv_cb_ptr != NULL) {
		fprintf(fp, "<span style='color:green'>%s</span>", func_deriv_cb_ptr->name);
		func_deriv_cb_ptr = func_deriv_cb_ptr->next;
	}
}

/* Writes a html row containing the adjoint variables into fp */
void adj_html_vars(FILE* fp, adj_adjointer* adjointer, int type)
{
  int i;
  char adj_name[ADJ_NAME_LEN];
  adj_variable adj_var;
  fprintf(fp, "<div style=\"height:150px\"></div>\n");
  fprintf(fp, "<tr>\n");
  for (i = 0; i < adjointer->nequations; i++)
  {
    adj_var = adjointer->equations[i].variable;
    adj_var.type = type;
    adj_variable_str(adj_var, adj_name, ADJ_NAME_LEN);
    fprintf(fp, "<th class=\"headercell box_rotate\"><div>%s</div></th>\n", adj_name);
  }
  fprintf(fp, "</tr>\n");
 }

/* Writes a html row containing the supplied equation into fp */
int adj_html_eqn(FILE* fp, adj_adjointer* adjointer, adj_equation adj_eqn, int diag_index, char* class)
{
  int i,k;
  char* row[adjointer->nequations];
  char* desc[adjointer->nequations];
  char buf[ADJ_NAME_LEN];
  int col, ierr;

  /* Allocate the strings for this row */
  for (i = 0; i < adjointer->nequations; ++i)
  {
    row[i] = malloc(ADJ_NAME_LEN*sizeof(char));
    ADJ_CHKMALLOC(row[i]);
    row[i][0]='\0';
    desc[i] = malloc(32*ADJ_NAME_LEN*sizeof(char));  // The description can become very long
    ADJ_CHKMALLOC(desc[i]);
    desc[i][0]='\0';
  }

  /* Fill in the data */
  for (i = 0; i < adj_eqn.nblocks; i++)
  {
    ierr = adj_html_find_column_index(adjointer, &adj_eqn.targets[i], &col);
    if (ierr != ADJ_OK)
      return ierr;

    strncpy(row[col], adj_eqn.blocks[i].name, 5);
    row[col][5]='\0';

    /* Fill in the description */
    strncpy(desc[col], adj_eqn.blocks[i].name, ADJ_NAME_LEN);

    strncat(desc[col], "\n\nTargets: ", ADJ_NAME_LEN);
    strncat(desc[col], adj_eqn.targets[i].name, ADJ_NAME_LEN);
    strncat(desc[col], "\nTimestep:", ADJ_NAME_LEN);
    snprintf(buf, ADJ_NAME_LEN, "%d", adj_eqn.targets[i].timestep);
    strncat(desc[col], buf, ADJ_NAME_LEN);
    strncat(desc[col], "\nIteration:", ADJ_NAME_LEN);
    snprintf(buf, ADJ_NAME_LEN, "%d", adj_eqn.targets[i].iteration);
    strncat(desc[col], buf, ADJ_NAME_LEN);
    strncat(desc[col], "\nCoefficient: ", ADJ_NAME_LEN);
    snprintf(buf, ADJ_NAME_LEN, "%f", adj_eqn.blocks[i].coefficient);
    strncat(desc[col], buf, ADJ_NAME_LEN);
    strncat(desc[col], "\nHermitian: ", ADJ_NAME_LEN);
    if (adj_eqn.blocks[i].hermitian==ADJ_TRUE)
      snprintf(buf, ADJ_NAME_LEN, "true");
    else
      snprintf(buf, ADJ_NAME_LEN, "false");
    strncat(desc[col], buf, ADJ_NAME_LEN);

    if (adj_eqn.blocks[i].has_nonlinear_block)
    {
      strncat(desc[col], "\nNonlinear Block: ", ADJ_NAME_LEN);
      strncat(desc[col], adj_eqn.blocks[i].nonlinear_block.name, ADJ_NAME_LEN);
      for (k=0; k<adj_eqn.blocks[i].nonlinear_block.ndepends; k++)
      {
        strncat(desc[col], " (Dependency: ", ADJ_NAME_LEN);
        strncat(desc[col], adj_eqn.blocks[i].nonlinear_block.depends[k].name, ADJ_NAME_LEN);
        strncat(desc[col], ":", ADJ_NAME_LEN);
        snprintf(buf, ADJ_NAME_LEN, "%d", adj_eqn.blocks[i].nonlinear_block.depends[k].timestep);
        strncat(desc[col], buf, ADJ_NAME_LEN);
        strncat(desc[col], ":", ADJ_NAME_LEN);
        snprintf(buf, ADJ_NAME_LEN, "%d", adj_eqn.blocks[i].nonlinear_block.depends[k].iteration);
        strncat(desc[col], buf, ADJ_NAME_LEN);
        strncat(desc[col], ")", ADJ_NAME_LEN);
      }

    }
  }
  /* Write it to file */
  adj_html_write_row(fp, row, desc, adjointer->nequations, diag_index, class);


  /* Tidy up */
  for (i = 0; i < adjointer->nequations; ++i)
  {
    free(row[i]);
    free(desc[i]);
  }
  return ADJ_OK;
}

/* Writes a html row containing the supplied adjoint equation into fp */
int adj_html_adjoint_eqn(FILE* fp, adj_adjointer* adjointer, adj_equation fwd_eqn, int diag_index, char* class)
{
  int i, j, k;
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
    ADJ_CHKMALLOC(row[i]);
    row[i][0]='\0';
    desc[i] = malloc(32*ADJ_NAME_LEN*sizeof(char)); // The description can become very long
    ADJ_CHKMALLOC(desc[i]);
    desc[i][0]='\0';
  }

  fwd_var = fwd_eqn.variable;
  ierr = adj_find_variable_data(&(adjointer->varhash), &fwd_var, &fwd_data);
  assert(ierr == ADJ_OK);

  /* Each targeting equation corresponds to one column in the considered row */
  for (i = 0; i < fwd_data->ntargeting_equations; i++)
  {
      adj_equation other_fwd_eqn;
      adj_variable other_adj_var;

      other_fwd_eqn = adjointer->equations[fwd_data->targeting_equations[i]];

      /* Find the index in other_fwd_eqn of the block fwd_var is multiplied with */
      for (j=0; j<other_fwd_eqn.nblocks; j++)
      {
        if (adj_variable_equal(&other_fwd_eqn.targets[j], &fwd_var, 1))
          break;
      }
      assert(j!=other_fwd_eqn.nblocks);

      /* find the column in which this blocks belongs */
      ierr = adj_html_find_column_index(adjointer, &other_fwd_eqn.variable, &col);
      if (ierr != ADJ_OK)
        return ierr;

      other_adj_var = other_fwd_eqn.targets[j];
      other_adj_var.type = ADJ_ADJOINT;

      /* Fill in the data */
      strncpy(row[col], other_fwd_eqn.blocks[j].name, 5);
      row[col][5]='\0';

      strncpy(desc[col], other_fwd_eqn.blocks[j].name, ADJ_NAME_LEN);
      strncat(desc[col], "\n\nTargets: ", ADJ_NAME_LEN);
      adj_variable_str(other_adj_var, buf, ADJ_NAME_LEN);
      strncat(desc[col], buf, ADJ_NAME_LEN);
      strncat(desc[col], "\nTimestep:", ADJ_NAME_LEN);
      snprintf(buf, ADJ_NAME_LEN, "%d", other_adj_var.timestep);
      strncat(desc[col], buf, ADJ_NAME_LEN);
      strncat(desc[col], "\nIteration:", ADJ_NAME_LEN);
      snprintf(buf, ADJ_NAME_LEN, "%d", other_adj_var.iteration);
      strncat(desc[col], buf, ADJ_NAME_LEN);
      strncat(desc[col], "\nCoefficient: ", ADJ_NAME_LEN);
      snprintf(buf, ADJ_NAME_LEN, "%f", other_fwd_eqn.blocks[j].coefficient);
      strncat(desc[col], buf, ADJ_NAME_LEN);
      strncat(desc[col], "\nHermitian: ", ADJ_NAME_LEN);
      /* We are printing the adjoint equation, therefore hermitian has to be NOT'ed */
      if (other_fwd_eqn.blocks[j].hermitian==ADJ_TRUE)
        snprintf(buf, ADJ_NAME_LEN, "false");
      else
        snprintf(buf, ADJ_NAME_LEN, "true");
      strncat(desc[col], buf, ADJ_NAME_LEN);
      if (other_fwd_eqn.blocks[j].has_nonlinear_block)
      {
        strncat(desc[col], "\nNonlinear Block: ", ADJ_NAME_LEN);
        strncat(desc[col], other_fwd_eqn.blocks[j].nonlinear_block.name, ADJ_NAME_LEN);
        for (k=0; k<other_fwd_eqn.blocks[j].nonlinear_block.ndepends; k++)
        {
          strncat(desc[col], " (Dependency: ", ADJ_NAME_LEN);
          strncat(desc[col], other_fwd_eqn.blocks[j].nonlinear_block.depends[k].name, ADJ_NAME_LEN);
          strncat(desc[col], ":", ADJ_NAME_LEN);
          snprintf(buf, ADJ_NAME_LEN, "%d", other_fwd_eqn.blocks[j].nonlinear_block.depends[k].timestep);
          strncat(desc[col], buf, ADJ_NAME_LEN);
          strncat(desc[col], ":", ADJ_NAME_LEN);
          snprintf(buf, ADJ_NAME_LEN, "%d", other_fwd_eqn.blocks[j].nonlinear_block.depends[k].iteration);
          strncat(desc[col], buf, ADJ_NAME_LEN);
          strncat(desc[col], ")", ADJ_NAME_LEN);
        }

      }

  }

  /* Write it to file */
  adj_html_write_row(fp, row, desc, adjointer->nequations, diag_index, class);

  /* Tidy up */
  for (i = 0; i < adjointer->nequations; ++i)
  {
    free(row[i]);
    free(desc[i]);
  }
  return ADJ_OK;
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

  adj_html_print_statistics(fp, adjointer);
  adj_html_print_callback_information(fp, adjointer);

  fprintf(fp, "<h1>Adjoint system</h1>\n");
  if (adjointer->nequations==0)
    return ADJ_OK;
  timestep = adjointer->equations[0].variable.timestep-1;
  iteration = adjointer->equations[0].variable.iteration - 1 ;

  adj_html_table_begin(fp, "");
  adj_html_vars(fp, adjointer, ADJ_ADJOINT);
  adj_html_table_end(fp);

  adj_html_table_begin(fp, "");
  for (i=0; i<adjointer->nequations; i++)
  {
    adj_eqn = adjointer->equations[i];
    if (timestep != adj_eqn.variable.timestep) {
    	ierr = adj_html_adjoint_eqn(fp, adjointer, adj_eqn, i, "new_timestep");
    }
    else {
      if (iteration != adj_eqn.variable.iteration)
    	  ierr = adj_html_adjoint_eqn(fp, adjointer, adj_eqn, i, "");
      else
    	  ierr = adj_html_adjoint_eqn(fp, adjointer, adj_eqn, i, "");
    }

    if (ierr != ADJ_OK)
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
  return ADJ_OK;
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

  adj_html_print_statistics(fp, adjointer);
  adj_html_print_callback_information(fp, adjointer);

  fprintf(fp, "<h1>Forward system</h1>\n");
  if (adjointer->nequations==0)
    return ADJ_OK;
  timestep = adjointer->equations[0].variable.timestep - 1;
  iteration = adjointer->equations[0].variable.iteration - 1 ;

  adj_html_table_begin(fp, "");
  adj_html_vars(fp, adjointer, ADJ_ADJOINT);
  adj_html_table_end(fp);

  adj_html_table_begin(fp, "");
  for (i=0; i<adjointer->nequations; i++)
  {
    adj_eqn = adjointer->equations[i];
    if (timestep != adj_eqn.variable.timestep) {
      ierr = adj_html_eqn(fp, adjointer, adj_eqn, i, "new_timestep");
    }
    else {
      if (iteration != adj_eqn.variable.iteration)
    	  ierr = adj_html_eqn(fp, adjointer, adj_eqn, i, "");
      else
    	  ierr = adj_html_eqn(fp, adjointer, adj_eqn, i, "");
    }


    if (ierr != ADJ_OK)
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
  return ADJ_OK;
}

int adj_adjointer_to_html(adj_adjointer* adjointer, char* filename, int type)
{
  if (type == ADJ_FORWARD)
    return adj_html_forward_system(adjointer, filename);
  else if(type == ADJ_ADJOINT)
    return adj_html_adjoint_system(adjointer, filename);
  else
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The type parameter must either be ADJ_FORWARD or ADJ_ADJOINT.");
    return ADJ_ERR_INVALID_INPUTS;
  }
}
