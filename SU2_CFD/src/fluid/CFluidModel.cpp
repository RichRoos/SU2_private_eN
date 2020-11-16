/*!
 * \file CFluidModel.cpp
 * \brief Source of the fluid model base class containing thermo-physical subroutines.
 * \author S.Vitale, M.Pini, G.Gori, A.Guardone, P.Colonna, T. Economon
 * \version 7.0.7 "Blackbird"
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

#include <utility>

#include "../../include/fluid/CFluidModel.hpp"
#include "../../include/fluid/CConstantConductivity.hpp"
#include "../../include/fluid/CConstantConductivityRANS.hpp"
#include "../../include/fluid/CConstantPrandtl.hpp"
#include "../../include/fluid/CConstantPrandtlRANS.hpp"
#include "../../include/fluid/CConstantViscosity.hpp"
#include "../../include/fluid/CPolynomialConductivity.hpp"
#include "../../include/fluid/CPolynomialConductivityRANS.hpp"
#include "../../include/fluid/CPolynomialViscosity.hpp"
#include "../../include/fluid/CSutherland.hpp"
#include "../../include/fluid/CConstantDiffusivity.hpp"
#include "../../include/fluid/CConstantSchmidtRANS.hpp"
#include "../../include/fluid/CConstantSchmidt.hpp"

void CFluidModel::SetLaminarViscosityModel(const CConfig* config) {
  switch (config->GetKind_ViscosityModel()) {
    case CONSTANT_VISCOSITY:
      LaminarViscosity = unique_ptr<CConstantViscosity>(new CConstantViscosity(config->GetMu_ConstantND()));
      break;
    case SUTHERLAND:
      LaminarViscosity = unique_ptr<CSutherland>(
          new CSutherland(config->GetMu_RefND(), config->GetMu_Temperature_RefND(), config->GetMu_SND()));
      break;
    case POLYNOMIAL_VISCOSITY:
      LaminarViscosity = unique_ptr<CPolynomialViscosity<N_POLY_COEFFS>>(
          new CPolynomialViscosity<N_POLY_COEFFS>(config->GetMu_PolyCoeffND()));
      break;
    case FLAMELET_VISC_MODEL:
      /* do nothing. Viscosity is obtained from the table and set in setTDState_T */
      break;
    default:
      SU2_MPI::Error("Viscosity model not available.", CURRENT_FUNCTION);
      break;
  }
}

void CFluidModel::SetThermalConductivityModel(const CConfig* config) {
  switch (config->GetKind_ConductivityModel()) {
    case CONSTANT_CONDUCTIVITY:
      if (config->GetKind_ConductivityModel_Turb() == CONSTANT_PRANDTL_TURB) {
        ThermalConductivity = unique_ptr<CConstantConductivityRANS>(
            new CConstantConductivityRANS(config->GetKt_ConstantND(), config->GetPrandtl_Turb()));
      } else {
        ThermalConductivity = unique_ptr<CConstantConductivity>(new CConstantConductivity(config->GetKt_ConstantND()));
      }
      break;
    case CONSTANT_PRANDTL:
      if (config->GetKind_ConductivityModel_Turb() == CONSTANT_PRANDTL_TURB) {
        ThermalConductivity = unique_ptr<CConstantPrandtlRANS>(
            new CConstantPrandtlRANS(config->GetPrandtl_Lam(), config->GetPrandtl_Turb()));
      } else {
        ThermalConductivity = unique_ptr<CConstantPrandtl>(new CConstantPrandtl(config->GetPrandtl_Lam()));
      }
      break;
    case POLYNOMIAL_CONDUCTIVITY:
      if (config->GetKind_ConductivityModel_Turb() == CONSTANT_PRANDTL_TURB) {
        ThermalConductivity = unique_ptr<CPolynomialConductivityRANS<N_POLY_COEFFS>>(
            new CPolynomialConductivityRANS<N_POLY_COEFFS>(config->GetKt_PolyCoeffND(), config->GetPrandtl_Turb()));
      } else {
        ThermalConductivity = unique_ptr<CPolynomialConductivity<N_POLY_COEFFS>>(
            new CPolynomialConductivity<N_POLY_COEFFS>(config->GetKt_PolyCoeffND()));
      }
      break;
    case FLAMELET_CONDUCT_MODEL:
      /* do nothing. Conductivity is obtained from the table and set in setTDState_T */
      break;
    default:
      SU2_MPI::Error("Conductivity model not available.", CURRENT_FUNCTION);
      break;
  }
}

void CFluidModel::SetMassDiffusivityModel (const CConfig* config) {
  switch (config->GetKind_DiffusivityModel()) {
    case CONSTANT_DIFFUSIVITY:
      MassDiffusivity = unique_ptr<CConstantDiffusivity>(new CConstantDiffusivity(config->GetDiffusivity_ConstantND()));
      break;
    case CONSTANT_SCHMIDT:
      if ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS)) {
        MassDiffusivity = unique_ptr<CConstantSchmidtRANS>(new CConstantSchmidtRANS(config->GetSchmidt_Lam(),config->GetSchmidt_Turb()));
      } else {
        MassDiffusivity = unique_ptr<CConstantSchmidt>(new CConstantSchmidt(config->GetSchmidt_Lam()));
      }
      break;
    case FLAMELET_DIFF_MODEL:
      /* do nothing. Diffusivity is obtained from the table and set in setTDState_T */
      break;
    default:
      SU2_MPI::Error("Diffusivity model not available.", CURRENT_FUNCTION);
      break;
  }
  
}