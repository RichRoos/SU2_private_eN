/*!
 * \file CScalarSolver.hpp
 * \brief Main subroutines for the  transported scalar model.
 * \author D. Mayer, T. Economon
 * \version 7.1.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "CSolver.hpp"
#include "../variables/CScalarVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"

/*!
 * \class CScalarSolver
 * \brief Main class for defining the solver for scalar transport eqns.
 * \ingroup Scalar_Model
 * \author T. Economon
 */
class CScalarSolver : public CSolver {
protected:
  su2double *FlowPrimVar_i,      /*!< \brief Store the flow solution at point i. */
  *FlowPrimVar_j,                /*!< \brief Store the flow solution at point j. */
  *lowerlimit,                   /*!< \brief contains lower limits for turbulence variables. */
  *upperlimit;                   /*!< \brief contains upper limits for turbulence variables. */
  vector<su2activematrix> Inlet_ScalarVars;  /*!< \brief scalar variables at inlet profiles */

 /* Sliding mesh variables */
  
  vector<su2matrix<su2double*> > SlidingState; // vector of matrix of pointers... inner dim alloc'd elsewhere (welcome, to the twilight zone)
  vector<vector<int> > SlidingStateNodes;

  unsigned long nMarker,         /*!< \brief Total number of markers using the grid information. */
  *nVertex;                      /*!< \brief Store nVertex at each marker for deallocation */
  unsigned short Inlet_Position; /*!< \brief Column index for scalar variables in inlet files. */
  su2double *Scalar_Inf;         /*!< \brief Array of far-field values for the scalar variables. */
  
  CScalarVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */
  
  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

private:

/*!
   * \brief Compute the viscous flux for the turbulent equation at a particular edge.
   * \param[in] iEdge - Edge for which we want to compute the flux
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Viscous_Residual(unsigned long iEdge,
                        CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *numerics,
                        CConfig *config);
  using CSolver::Viscous_Residual; /*--- Silence warning ---*/

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CScalarSolver(void);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CScalarSolver(void);
  
  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CScalarSolver(CGeometry* geometry, CConfig *config);
  
  /*!
   * \brief Compute the spatial integration using a upwind scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  
  void Upwind_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics **numerics_container,
                       CConfig *config,
                       unsigned short iMesh) override;
  
  
  
  /*!
   * \brief Prepare an implicit iteration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void PrepareImplicitIteration(CGeometry *geometry, CSolver** solver_container, CConfig *config) final;

  /*!
   * \brief Complete an implicit iteration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void CompleteImplicitIteration(CGeometry *geometry, CSolver** solver_container, CConfig *config) final;


  /*!
   * \brief Update the solution using an implicit scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitEuler_Iteration(CGeometry  *geometry,
                               CSolver   **solver_container,
                               CConfig    *config) override;
  
  /*!
   * \brief Impose the Symmetry Plane boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Sym_Plane(CGeometry *geometry,
                    CSolver **solver_container,
                    CNumerics *conv_numerics,
                    CNumerics *visc_numerics,
                    CConfig *config,
                    unsigned short val_marker) override;
  
  /*!
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Euler_Wall(CGeometry *geometry,
                     CSolver **solver_container,
                     CNumerics *conv_numerics,
                     CNumerics *visc_numerics,
                     CConfig *config,
                     unsigned short val_marker) override;
  
  /*!
   * \brief Impose the Navier-Stokes wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container,
                        CNumerics *conv_numerics, CNumerics *visc_numerics,
                        CConfig *config, unsigned short val_marker) override;
    
  /*!
   * \brief Impose the Far Field boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Far_Field(CGeometry       *geometry, 
                    CSolver        **solver_container,
                    CNumerics       *conv_numerics,
                    CNumerics       *visc_numerics,
                    CConfig         *config, 
                    unsigned short   val_marker) override;

  /*!
   * \brief Impose a periodic boundary condition by summing contributions from the complete control volume.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                   CNumerics *numerics, CConfig *config) override;



  /*!
   * \brief Set the total residual adding the term that comes from the Dual Time-Stepping Strategy.
   * \param[in] geometry - Geometric definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void SetResidual_DualTime(CGeometry       *geometry,
                            CSolver        **solver_container,
                            CConfig         *config,
                            unsigned short   iRKStep,
                            unsigned short   iMesh,
                            unsigned short   RunTime_EqSystem) override;
  
  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry  **geometry,
                   CSolver   ***solver,
                   CConfig     *config,
                   int          val_iter,
                   bool         val_update_geo) override;
  
  /*!
   * \brief Set custom scalar variables at the vertex of an inlet.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] iDim - Index of the scalar variable
   * \param[in] val_turb_var - Value of the turbulence variable to be used.
   */
  void SetInlet_ScalarVar(unsigned short val_marker,
                          unsigned long val_vertex,
                          unsigned short val_dim,
                          su2double val_scalar_var);
  
  /*!
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  inline void SetFreeStream_Solution(const CConfig *config) override {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
      for (unsigned short iVar = 0; iVar < nVar; iVar++)
        nodes->SetSolution(iPoint, iVar, Scalar_Inf[iVar]);
  }
  
  /*!
   * \brief Store of a set of provided inlet profile values at a vertex.
   * \param[in] val_inlet - vector containing the inlet values for the current vertex.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   */
  void SetInletAtVertex(const su2double *val_inlet, 
                        unsigned short  iMarker,
                        unsigned long   iVertex) override;
  
  /*!
   * \brief Get the set of values imposed at an inlet.
   * \param[in] val_inlet - vector returning the inlet values for the current vertex.
   * \param[in] val_inlet_point - Node index where the inlet is being set.
   * \param[in] val_kind_marker - Enumerated type for the particular inlet type.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param config - Definition of the particular problem.
   * \return Value of the face area at the vertex.
   */
  su2double GetInletAtVertex(su2double       *val_inlet,
                             unsigned long    val_inlet_point,
                             unsigned short   val_kind_marker,
                             string           val_marker,
                             const CGeometry *geometry,
                             const CConfig   *config) const override;
  
  /*!
   * \brief Set a uniform inlet profile
   *
   * The values at the inlet are set to match the values specified for
   * inlets in the configuration file.
   *
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   */
  void SetUniformInlet(const CConfig *config, unsigned short iMarker) override;

  /*!
   * \brief Get the value of the scalar variables at the far-field.
   * \return Value of the scalar variables at the far-field.
   */
  inline su2double GetScalar_Inf(unsigned short val_ivar) { return Scalar_Inf[val_ivar]; }
  
};