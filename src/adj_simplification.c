#include "libadjoint/adj_simplification.h"

int adj_simplify_derivatives(adj_adjointer* adjointer, int ninput, adj_nonlinear_block_derivative* input, int* noutput, adj_nonlinear_block_derivative** output)
{
  /* During adj_get_adjoint_equation, we construct a list of derivatives to be computed; these derivatives
     are stored in an adj_nonlinear_block_derivative. They take the form

              ([ dV ]   )*
     factor * ([ -- ] c )  lambda
              ([ du ]   )

     Now, if two adj_nonlinear_block_derivatives share the same u and the same c, they can be merged, meaning
     we only have to call the (potentially expensive) derivative computation routine.

     Doing this in a naive way requires you to compare every derivative with every other derivative.
     That's quadratic, and bad; however, it's much easier to code, so that's what we're going to do
     for now. Typically we will only have ~4-5 of these to process at once anyhow, so it probably
     won't ever be a problem.

     For your edification, I will record the linear algorithm here, so that if it ever becomes a problem,
     you'll know what to do!

     The linear solution involves creating a hash table that maps
     (variable to differentiate, operator name, operator dependencies, context)
     to
     (a linked list of adj_nonlinear_block_derivatives that all have the same properties).

     Use uthash, in a similar way to the variable hash table that forms the heart of the
     adjointer.

     Once you've finished hashing, walk the lists.
     For each list:
       merge all the entries in the list.

     And you're done. If it's your job to code it, have fun! -- pef */

  int i;
  int j;
  int count;
  int merged[ninput];
  int discarded[ninput];
  int return_ierr;
  adj_nonlinear_block_derivative* copy;

  memset(merged, ADJ_FALSE, ninput * sizeof(int));
  memset(discarded, ADJ_FALSE, ninput * sizeof(int));

  copy = (adj_nonlinear_block_derivative*) malloc(ninput * sizeof(adj_nonlinear_block_derivative));
  ADJ_CHKMALLOC(copy);
  memcpy(copy, input, ninput * sizeof(adj_nonlinear_block_derivative));

  for (i = 0; i < ninput; i++)
  {
    for (j = i+1; j < ninput; j++)
    {
      if (discarded[j]) continue;
      if (adj_simplification_compare(copy[i], copy[j]))
      {
        adj_simplification_merge(adjointer, &copy[i], &copy[j], merged[i]);
        discarded[j] = ADJ_TRUE;
        merged[i] = ADJ_TRUE;
      }
    }
  }

  count = 0;
  for (i = 0; i < ninput; i++)
    if (!discarded[i]) count++;

  *noutput = count;
  *output = (adj_nonlinear_block_derivative*) malloc(count * sizeof(adj_nonlinear_block_derivative));
  ADJ_CHKMALLOC(*output);

  /* Compact here */
  j = 0;
  for (i = 0; i < ninput; i++)
  {
    if (!discarded[i])
    {
      if (merged[i])
      {
        (*output)[j] = copy[i];
      }
      else
      {
        /* We need to make a fresh copy, because adj_get_adjoint_equation will be deallocating these
           very soon. */
        int ierr;
        ierr = adj_create_nonlinear_block_derivative(adjointer, copy[i].nonlinear_block, (adj_scalar) 1.0, copy[i].variable, copy[i].contraction, copy[i].hermitian, &(*output)[j]);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
      }
      j++;
    }
  }
  free(copy);

  if (ninput > 45) /* we did more than 1000 comparisons */
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "adj_simplify_derivatives contains a known quadratic loop.\nA linear algorithm is known and discussed in the comments.\nYou may wish to consider implementing it, since your problem\nis getting quite large (45 blocks in a row? really?).");
    return_ierr = ADJ_WARN_NOT_IMPLEMENTED;
  }
  else
  {
    return_ierr = ADJ_OK;
  }

  return return_ierr;
}

int adj_simplification_compare(adj_nonlinear_block_derivative d1, adj_nonlinear_block_derivative d2)
{
  return (adj_variable_equal(&d1.variable, &d2.variable, 1)) &&
         (d1.nonlinear_block.ndepends == d2.nonlinear_block.ndepends) &&
         (adj_variable_equal(d1.nonlinear_block.depends, d2.nonlinear_block.depends, d1.nonlinear_block.ndepends)) &&
         (strncmp(d1.nonlinear_block.name, d2.nonlinear_block.name, ADJ_NAME_LEN) == 0) &&
         (d1.nonlinear_block.context == d2.nonlinear_block.context) &&
         (d1.hermitian == d2.hermitian);
}

void adj_simplification_merge(adj_adjointer* adjointer, adj_nonlinear_block_derivative* d1, adj_nonlinear_block_derivative* d2, int merged)
{
  adj_vector new_contraction;

  adjointer->callbacks.vec_duplicate(d1->contraction, &new_contraction);
  adjointer->callbacks.vec_axpy(&new_contraction, d1->nonlinear_block.coefficient, d1->contraction);
  adjointer->callbacks.vec_axpy(&new_contraction, d2->nonlinear_block.coefficient, d2->contraction);
  d1->nonlinear_block.coefficient = (adj_scalar) 1.0;

  if (merged) /* We created this earlier as a result of a merge */
  {
    adjointer->callbacks.vec_destroy(&d1->contraction);
  }
  d1->contraction = new_contraction;
}
