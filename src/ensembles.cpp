/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/





#include <ensembles.hpp>
#include <XForm.hpp>
#include <AtomicGroup.hpp>
#include <Trajectory.hpp>
#include <alignment.hpp>

namespace loos {

  // Assume all groups are already sorted or matched...

  AtomicGroup averageStructure(const std::vector<AtomicGroup>& ensemble) {
    AtomicGroup avg = ensemble[0].copy();

    // First, zap our coords...
    int n = avg.size();
    int i;
    for (i=0; i<n; i++)
      avg[i]->coords() = GCoord(0.0, 0.0, 0.0);

    // Now, accumulate...
    std::vector<AtomicGroup>::const_iterator j;
    for (j = ensemble.begin(); j != ensemble.end(); j++) {
      for (i = 0; i<n; i++)
        avg[i]->coords() += (*j)[i]->coords();
    }

    for (i=0; i<n; i++)
      avg[i]->coords() /= ensemble.size();

    avg.removePeriodicBox();
    return(avg);
  }



  AtomicGroup averageStructure(const std::vector<AtomicGroup>& ensemble, const std::vector<XForm>& xforms) {
    if (xforms.size() != ensemble.size())
      throw(LOOSError("Transforms do not match the passed ensemble in loos::averageStructure()"));

    AtomicGroup avg = ensemble[0].copy();

    // First, zap our coords...
    int n = avg.size();
    int i;
    for (i=0; i<n; i++)
      avg[i]->coords() = GCoord(0.0, 0.0, 0.0);

    // Now, accumulate...
    for (uint j=0; j<ensemble.size(); ++j) {
      AtomicGroup frame = ensemble[j].copy();
      frame.applyTransform(xforms[j]);
      for (i = 0; i<n; i++)
        avg[i]->coords() += frame[i]->coords();
    }
    
    for (i=0; i<n; i++)
      avg[i]->coords() /= ensemble.size();

    avg.removePeriodicBox();
    return(avg);
  }




  AtomicGroup averageStructure(const AtomicGroup& g, const std::vector<XForm>& xforms, pTraj& traj, const std::vector<uint>& frame_indices) {
    AtomicGroup avg = g.copy();
    AtomicGroup frame = g.copy();
    int n = avg.size();
    for (int i=0; i<n; i++)
      avg[i]->coords() = GCoord(0.0, 0.0, 0.0);
    
    uint tn = traj->nframes();
    uint fn = frame_indices.size();
    
    // Safety check...
    if (fn != xforms.size())
      throw(std::runtime_error("Mismatch in number of frames in the trajectory requested and passed transforms for loos::averageStructure()"));
    
    for (uint j=0; j<fn; ++j) {
      if (frame_indices[j] >= tn)
        throw(std::runtime_error("Frame index exceeds trajectory size"));
      
      traj->readFrame(frame_indices[j]);
      traj->updateGroupCoords(frame);
      frame.applyTransform(xforms[j]);
      for (int i=0; i<n; i++)
        avg[i]->coords() += frame[i]->coords();
    }
    
    for (int i=0; i<n; i++)
      avg[i]->coords() /= fn;

    avg.removePeriodicBox();
    return(avg);
  }
  
  
  
  AtomicGroup averageStructure(const AtomicGroup& g, const std::vector<XForm>& xforms, pTraj& traj) {
    uint nx = xforms.size();
    uint nt = traj->nframes();
    
    if (nt != nx)
      throw(std::runtime_error("Mismatch in number of frames in the trajectory and passed transforms"));
    std::vector<uint> frame_indices;
    for (uint i=0; i<nt; ++i)
      frame_indices.push_back(i);

    return(averageStructure(g, xforms, traj, frame_indices));
  }



  void applyTransforms(std::vector<AtomicGroup>& ensemble, std::vector<XForm>& xforms) {
    uint n = ensemble.size();
    if (n != xforms.size())
      throw(std::runtime_error("Mismatch in the size of the ensemble and the transformations"));

    for (uint i=0; i<n; ++i)
      ensemble[i].applyTransform(xforms[i]);
  }



  void readTrajectory(std::vector<AtomicGroup>& ensemble, const AtomicGroup& model, pTraj trajectory) {
    AtomicGroup clone = model.copy();
    
    while (trajectory->readFrame()) {
      trajectory->updateGroupCoords(clone);
      AtomicGroup frame = clone.copy();
      ensemble.push_back(frame);
    }
  }


  void readTrajectory(std::vector<AtomicGroup>& ensemble, const AtomicGroup& model, pTraj trajectory, std::vector<uint>& frames) {
    AtomicGroup clone = model.copy();
    
    std::vector<uint>::iterator i;
    for (i = frames.begin(); i != frames.end(); ++i) {
      if (*i >= trajectory->nframes())
        throw(std::runtime_error("Frame index exceeds trajectory size in readTrajectory()"));
      trajectory->readFrame(*i);
      trajectory->updateGroupCoords(clone);
      AtomicGroup frame = clone.copy();
      ensemble.push_back(frame);
    }
  }



  RealMatrix extractCoords(const std::vector<AtomicGroup>& ensemble) {
    uint n = ensemble.size();
    uint m = ensemble[0].size();
    RealMatrix M(3*m, n);

    for (uint i=0; i<n; ++i)
      for (uint j=0; j<m; ++j) {
        M(3*j, i) = ensemble[i][j]->coords().x();
        M(3*j+1, i) = ensemble[i][j]->coords().y();
        M(3*j+2, i) = ensemble[i][j]->coords().z();
      }

    return(M);
  }


  RealMatrix extractCoords(const std::vector<AtomicGroup>& ensemble, const std::vector<XForm>& xforms) {
    uint n = ensemble.size();

    if (n != xforms.size())
      throw(std::runtime_error("Mismatch between the size of the ensemble and the transformations"));

    uint m = ensemble[0].size();
    RealMatrix M(3*m, n);

    for (uint i=0; i<n; ++i) {
      GMatrix W = xforms[i].current();

      for (uint j=0; j<m; ++j) {
        GCoord c = W * ensemble[i][j]->coords();
        M(3*j, i) = c.x();
        M(3*j+1, i) = c.y();
        M(3*j+2, i) = c.z();
      }
    }

    return(M);
  }


  void subtractAverage(RealMatrix& M) {
    uint m = M.rows();
    uint n = M.cols();
    double *avg = new double[m];
    
    for (uint j=0; j<m; ++j)
      avg[j] = 0.0;

    for (uint i=0; i<n; ++i)
      for (uint j=0; j<m; ++j)
        avg[j] += M(j, i);

    for (uint j=0; j<m; ++j)
      avg[j] /= n;

    for (uint i=0; i<n; ++i)
      for (uint j=0; j<m; ++j)
        M(j,i) -= avg[j];
  }


  boost::tuple<RealMatrix, RealMatrix, RealMatrix> svd(std::vector<AtomicGroup>& ensemble, bool align) {

    RealMatrix M;
    if (align)
      iterativeAlignment(ensemble);
    
    M = extractCoords(ensemble);

    subtractAverage(M);
    boost::tuple<RealMatrix, RealMatrix, RealMatrix> res = Math::svd(M);
    return(res);
  }


  void appendCoords(std::vector< std::vector<double> >& M, AtomicGroup& model, pTraj& traj, const std::vector<uint>& indices, const bool updates = false) {
    
    uint l = indices.size();
    uint n = model.size();
    uint offset = M.size();

    M.resize(offset + l, std::vector<double>(3*n));
    
    PercentProgressWithTime watcher;
    PercentTrigger trigger(0.1);
    ProgressCounter<PercentTrigger, EstimatingCounter> slayer(trigger, EstimatingCounter(l));
    
    if (updates) {
      slayer.attach(&watcher);
      slayer.start();
    }
  
    for (uint j=0; j<l; ++j) {
      traj->readFrame(indices[j]);
      traj->updateGroupCoords(model);
      if (updates)
        slayer.update();
      for (uint i=0; i<n; ++i) {
        GCoord c = model[i]->coords();
        M[j + offset][i*3] = c.x();
        M[j + offset][i*3+1] = c.y();
        M[j + offset][i*3+2] = c.z();
      }
    }
    
    if (updates)
      slayer.finish();

  }

  


  std::vector< std::vector<double> > readCoords(AtomicGroup& model, pTraj& traj, const std::vector<uint>& indices, const bool updates = false) {
    
    std::vector< std::vector<double> > M;
    appendCoords(M, model, traj, indices, updates);
    return M;
  }




  
}
