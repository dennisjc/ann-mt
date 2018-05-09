#include <iostream>
#include <fstream>
#include <string>
using namespace std;
int main (){
  // test reading from a file (inspired by:
  //http://stackoverflow.com/questions/14516915/read-numeric-data-from-a-text-file-in-c

  string firstline;
  string d;
  
  ifstream myfile;
  myfile.open ("test_model.00");
  getline(myfile,firstline);
  cout << firstline << endl;
  for(int i=0; i<2;i++){
    myfile>>d;
    cout<<d<<" ";
  }
  getline(myfile,firstline);
  cout << firstline;
  myfile>>d;
  cout << d;
  myfile.close();
  return 0;

}
