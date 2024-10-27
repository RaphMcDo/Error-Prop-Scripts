
using namespace Eigen;

template <class Type>
matrix <Type> DistMat2d(matrix <Type> locations) {
  
  int n_obs = locations.col(0).size();
  
  matrix <Type> distance(n_obs,n_obs);
  
  for (int i=0; i < n_obs; i++){
    distance(i,i) = Type(0.0);
    for (int j=i+1; j < n_obs; j++) {
      distance(i,j) = sqrt(pow(locations(i,0)-locations(j,0),2)+pow(locations(i,1)-locations(j,1),2));
      distance(j,i) = distance(i,j);
    }
  }
  
  return distance;
  
}


template <class Type>
matrix <Type> DistMat2d_aniso(matrix <Type> locations, matrix <Type> H_aniso) {
  
  int n_obs = locations.col(0).size();
  
  matrix <Type> distance(n_obs,n_obs);
  
  matrix <Type> mod_loc(n_obs,n_obs);
  
  for (int i=0; i < n_obs; i++){
    matrix <Type> temp_mat1(2,1);
    temp_mat1(0,0) = locations(i,0);
    temp_mat1(1,0) = locations(i,1);
    matrix <Type> temp_mat2(2,1); 
    temp_mat2 = H_aniso * temp_mat1;
    mod_loc(i,0) = temp_mat2(0,0);
    mod_loc(i,1) = temp_mat2(1,0);
  }
  
  for (int i=0; i < n_obs; i++){
    distance(i,i) = Type(0.0);
    for (int j=i+1; j < n_obs; j++) {
      distance(i,j) = sqrt(pow(mod_loc(i,0)-mod_loc(j,0),2)+pow(mod_loc(i,1)-mod_loc(j,1),2));
      distance(j,i) = distance(i,j);
    }
  }
  
  return distance;
  
}



