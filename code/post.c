#include "header.h"

int main(int argc, char **argv)
{

  post_str      post;
  geom_str      geom;

  /* --------------------------------------------- */

  MPI_Init(&argc, &argv);

  post_param(&argc, &argv, &post, &geom);
  //printf("OK\n");   MPI_Finalize(); exit(0); 


  timers_init(post.runname, 0);


  if( strcmp(post.type,"blur") == 0 ) post_blur(&post, &geom);
  else if ( strcmp(post.type,"bumps") == 0 ) post_bumps(&post, &geom);
  else if ( strcmp(post.type,"corrfun") == 0 ) post_corrfun(&post, &geom);
  else if ( strcmp(post.type,"corrfunx") == 0 ) post_corrfunx(&post, &geom);
  else if ( strcmp(post.type,"condensate") == 0 ) post_condensate(&post, &geom);
  else printf("unknown post-processing type:  %s\n", post.type);

  MPI_Finalize();
  exit(0); 
}

