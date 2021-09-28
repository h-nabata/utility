int digitcount(int number){
  int di = 0;
  int tempval = number;
  while( tempval != 0 ){
    tempval = tempval / 10;
    di++;
  }
  return di;
}
