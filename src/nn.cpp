#include "Progress.h"
#include <Rcpp.h>
#include "ANN.h"     // ANN library header

using namespace Rcpp;

NumericVector Cquantile(NumericVector x, double q) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y[x.size()*(q - 0.000000001)];
}

// [[Rcpp::export]]
NumericVector C_lasmerge_neighbors(S4 las, S4 las_from, const unsigned int k, const std::string func, const std::string var)
{
  DataFrame data = as<DataFrame>(las.slot("data"));
  DataFrame data_from = as<DataFrame>(las_from.slot("data"));

  const int d=3; //dimension of the point to index
  const int nd=data_from.nrow();
  const int nq=data.nrow();

  NumericVector X = data_from["X"];
  NumericVector Y = data_from["Y"];
  NumericVector Z = data_from["Z"];
  NumericVector X_query = data["X"];
  NumericVector Y_query = data["Y"];
  NumericVector Z_query = data["Z"];

  ANNkd_tree	*the_tree;	// Search structure

  ANNpointArray data_pts 	= annAllocPts(nd,d);		// Allocate data points
  ANNidxArray nn_idx 		= new ANNidx[k];		// Allocate near neigh indices
  ANNdistArray dists 		= new ANNdist[k];		// Allocate near neighbor dists

  Progress pb(nd, "Building index: ");
  // now construct the points
  for(int i = 0; i < nd; i++)
  {
    pb.check_abort();
    pb.increment();
    data_pts[i][0]=X[i];
    data_pts[i][1]=Y[i];
    data_pts[i][2]=Z[i];
  }
  the_tree = new ANNkd_tree( data_pts, nd, d);

  // return values here

  //now iterate over query points
  ANNpoint pq = annAllocPt(d);
  NumericVector V = data_from[var];

  Progress pb2(nq, "Computing neighbors values: ");

  NumericVector output(nq);
  if (k == 1) {
    for(int i = 0; i < nq; i++)	// Run all query points against tree
    {
      pb2.check_abort();
      pb2.increment();

      // read coords of current query point
      pq[0]=X_query(i);
      pq[1]=Y_query(i);
      pq[2]=Z_query(i);

      the_tree->annkSearch(	// search
          pq,	// query point
          k,		// number of near neighbors
          nn_idx,		// nearest neighbors (returned)
          dists,		// distance (returned)
          0.0);	// error bound

      output[i] = V[nn_idx[0]];
    }
  } else {
    for(int i = 0; i < nq; i++)	// Run all query points against tree
    {
      pb2.check_abort();
      pb2.increment();

      // read coords of current query point
      pq[0]=X_query(i);
      pq[1]=Y_query(i);
      pq[2]=Z_query(i);

      the_tree->annkSearch(	// search
          pq,	// query point
          k,		// number of near neighbors
          nn_idx,		// nearest neighbors (returned)
          dists,		// distance (returned)
          0.0);	// error bound

      NumericVector temp_neighbors(k);
      for (unsigned int j = 0 ; j < k ; j++)
      {
        temp_neighbors[j] = V[nn_idx[j]];
      }

      if (func == "var")
        output[i] = Rcpp::var(temp_neighbors);
      else if (func == "sd")
        output[i] = sd(temp_neighbors);
      else if (func == "median")
        output[i] = median(temp_neighbors);
      else if (func == "min")
        output[i] = min(temp_neighbors);
      else if (func == "max")
        output[i] = max(temp_neighbors);
      else if (func == "q1")
        output[i] = Cquantile(temp_neighbors,0.25)[1];
      else if (func == "q2")
        output[i] = Cquantile(temp_neighbors,0.5)[1];
      else if (func == "q3")
        output[i] = Cquantile(temp_neighbors,0.75)[1];
      else
        output[i] = mean(temp_neighbors);
    }
  }

  annDeallocPt(pq);
  annDeallocPts(data_pts);
  delete the_tree;
  delete [] nn_idx;
  delete [] dists;

  return (output);
}


// [[Rcpp::export]]
NumericVector C_lasmerge_neighbors_Rfunc(S4 las, S4 las_from, const unsigned int k, Rcpp::Function func, const std::string var)
{
  DataFrame data = as<DataFrame>(las.slot("data"));
  DataFrame data_from = as<DataFrame>(las_from.slot("data"));

  const int d=3; //dimension of the point to index
  const int nd=data_from.nrow();
  const int nq=data.nrow();

  NumericVector X = data_from["X"];
  NumericVector Y = data_from["Y"];
  NumericVector Z = data_from["Z"];
  NumericVector X_query = data["X"];
  NumericVector Y_query = data["Y"];
  NumericVector Z_query = data["Z"];

  ANNkd_tree	*the_tree;	// Search structure

  ANNpointArray data_pts 	= annAllocPts(nd,d);		// Allocate data points
  ANNidxArray nn_idx 		= new ANNidx[k];		// Allocate near neigh indices
  ANNdistArray dists 		= new ANNdist[k];		// Allocate near neighbor dists

  Progress pb(nd, "Building index: ");
  // now construct the points
  for(int i = 0; i < nd; i++)
  {
    pb.check_abort();
    pb.increment();
    data_pts[i][0]=X[i];
    data_pts[i][1]=Y[i];
    data_pts[i][2]=Z[i];
  }
  the_tree = new ANNkd_tree( data_pts, nd, d);

  // return values here

  //now iterate over query points
  ANNpoint pq = annAllocPt(d);
  NumericVector V = data_from[var];

  Progress pb2(nq, "Computing neighbors values: ");

  NumericVector output(nq);
  if (k == 1) {
    for(int i = 0; i < nq; i++)	// Run all query points against tree
    {
      pb2.check_abort();
      pb2.increment();

      // read coords of current query point
      pq[0]=X_query(i);
      pq[1]=Y_query(i);
      pq[2]=Z_query(i);

      the_tree->annkSearch(	// search
          pq,	// query point
          k,		// number of near neighbors
          nn_idx,		// nearest neighbors (returned)
          dists,		// distance (returned)
          0.0);	// error bound

      output[i] = V[nn_idx[0]];
    }
  } else {
    for(int i = 0; i < nq; i++)	// Run all query points against tree
    {
      pb2.check_abort();
      pb2.increment();

      // read coords of current query point
      pq[0]=X_query(i);
      pq[1]=Y_query(i);
      pq[2]=Z_query(i);

      the_tree->annkSearch(	// search
          pq,	// query point
          k,		// number of near neighbors
          nn_idx,		// nearest neighbors (returned)
          dists,		// distance (returned)
          0.0);	// error bound

      NumericVector temp_neighbors(k);
      for (unsigned int j = 0 ; j < k ; j++)
      {
        temp_neighbors[j] = V[nn_idx[j]];
      }
      NumericVector result = func(temp_neighbors);
      output[i] = result[0];
    }
  }

  annDeallocPt(pq);
  annDeallocPts(data_pts);
  delete the_tree;
  delete [] nn_idx;
  delete [] dists;

  return (output);
}
