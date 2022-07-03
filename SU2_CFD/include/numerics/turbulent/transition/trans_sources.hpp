/*!
 * \file trans_sources.hpp
 * \brief Numerics classes for integration of source terms in turbulence problems.
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../scalar/scalar_sources.hpp"


/*!
 * \class CSourcePieceWise_TranEN
 * \brief Class for integrating the source terms of the Menter SST turbulence model equations.
 * \ingroup SourceDiscr
 * \author S. Kang.
 */
template <class FlowIndices>
class CSourcePieceWise_TransEN final : public CNumerics {
 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */

  su2double g_eff_i,
  g_eff_j,
  g_sep_i,
  g_sep_j;
  /*--- eN Closure constants ---*/

  su2double Vorticity;
  su2double Residual[2];
  su2double* Jacobian_i[2];
  su2double Jacobian_Buffer[4]; /// Static storage for the Jacobian (which needs to be pointer for return type).

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] constants - SST model constants.
   * \param[in] val_kine_Inf - Freestream k, for SST with sustaining terms.
   * \param[in] val_omega_Inf - Freestream w, for SST with sustaining terms.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TransEN(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) 
      : CNumerics(val_nDim, 2, config),
        idx(val_nDim, config->GetnSpecies()) {

    
    /*--- "Allocate" the Jacobian using the static buffer. ---*/
    Jacobian_i[0] = Jacobian_Buffer;
    Jacobian_i[1] = Jacobian_Buffer + 2;
  }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override {
  AD::StartPreacc();
  AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(ScalarVar_i, nVar);
  AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(TransVar_i, nVar);
  AD::SetPreaccIn(TransVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);
  AD::SetPreaccIn(&V_i[idx.Velocity()], nDim);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+idx.Velocity(), nDim);
  AD::SetPreaccIn(Vorticity_i, 3);

  unsigned short iDim;
  /*
  su2double Corr_Rec, Corr_F_length, Corr_Ret, F_length;
  su2double F_onset1, F_onset2, F_onset3, F_onset;
  su2double f_turb, Re_v, f_sub, r_omega;
  su2double Pg, Dg, Pthetat;
  su2double diverg;
  su2double lambda, Tu;
  su2double f_lambda, Corr_Ret_lim, R_t;
  su2double vel_u, vel_v, vel_w, Velocity_Mag = 0.0, du_ds;
  su2double theta,  time_scale, f_theta;  
  su2double dU_dx, dU_dy, dU_dz = 0.0;
  su2double theta_bl, delta_bl, delta;
  su2double f_reattach, f_wake, re_omega, var1, var2;
  */

  su2double VorticityMag = sqrt(Vorticity_i[0]*Vorticity_i[0] +
                                Vorticity_i[1]*Vorticity_i[1] +
                                Vorticity_i[2]*Vorticity_i[2]);
  
  diverg = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) diverg += PrimVar_Grad_i[iDim+idx.Velocity()][iDim];

  vel_u = 0.0, vel_v = 0.0, vel_w = 0.0;
  if(nDim ==2){
      vel_u = V_i[idx.Velocity()];
      vel_v = V_i[1+idx.Velocity()];
  }
  else if(nDim ==3){
      vel_u = V_i[idx.Velocity()];
      vel_v = V_i[1+idx.Velocity()];
      vel_w = V_i[2+idx.Velocity()];
  }

  Velocity_Mag = sqrt(vel_u*vel_u + vel_v*vel_v + vel_w*vel_w);

  AD::SetPreaccIn(V_i[idx.Density()], V_i[idx.LaminarViscosity()], V_i[idx.EddyViscosity()]);

  Density_i = V_i[idx.Density()];
  Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];
  Eddy_Viscosity_i = V_i[idx.EddyViscosity()];

  Tu = 100.0*sqrt( 2.0 * ScalarVar_i[0] / 3.0 ) / Velocity_Mag;

  Residual[0] = 0.0;       Residual[1] = 0.0;
  Jacobian_i[0][0] = 0.0;  Jacobian_i[0][1] = 0.0;
  Jacobian_i[1][0] = 0.0;  Jacobian_i[1][1] = 0.0;
  
  if (dist_i > 1e-10) {

  }
  
  AD::SetPreaccOut(intermittency_sep_i);
  AD::SetPreaccOut(Residual, nVar);
  AD::EndPreacc();

  return ResidualType<>(Residual, Jacobian_i, nullptr);
  }
};

