int strcount(char *str1, char *str2) {
  char *start_ptr;  /* search begin point */
  char *search_ptr; /* search end point */
  int occurrences_number = 0;

  start_ptr = str1;
  search_ptr = strstr(start_ptr, str2);
  while( search_ptr != NULL ){
    start_ptr = search_ptr + strlen(str2);
    occurrences_number++;
    search_ptr = strstr(start_ptr, str2);
  }
  return occurrences_number;
}
