#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"
//#define MIN(a,b) (((a)<(b))?(a):(b))
//#define MAX(a,b) (((a)>(b))?(a):(b))


extern "C" void ApplyFakeTmunu(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // number of grid points
  const int nx = cctk_lsh[0];
  const int ny = cctk_lsh[1];
  const int nz = cctk_lsh[2];
  const double time = cctk_time; 
 
#pragma omp parallel for
    for (int k=0;k<nz;k++) {
      for (int j=0;j<ny;j++) {
	      for (int i=0;i<nx;i++) {

          // indices
          const int i3D = CCTK_GFINDEX3D(cctkGH,i,j,k);
          const int xi3D = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
          const int yi3D = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1);
          const int zi3D = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2);

          CCTK_REAL phidot, phi, rtraj;
          /* Select the trajectory */
	  //if (FakeTmunu_traj == "circular orbit") { 
          //if ( CCTK_EQUALS(FakeTmunu_traj, "circular orbit") ) {
	  //      phidot = FakeTmunu_omega;
	  //      phi = FakeTmunu_omega * time;
          //      rtraj = FakeTmunu_r;
          //}
          phidot = FakeTmunu_omega;
          phi = FakeTmunu_omega * time;
          rtraj = FakeTmunu_r;


          CCTK_REAL xtraj = rtraj * cos(phi);
	  CCTK_REAL ytraj = rtraj * sin(phi);
	  CCTK_REAL ztraj = 0.0;

	  CCTK_REAL vconx_traj =  rtraj * phidot * sin(phi);
	  CCTK_REAL vcony_traj = -rtraj * phidot * cos(phi);
	  CCTK_REAL vconz_traj =  0.0;

	  CCTK_REAL vcovx_traj =  gxx[i3D]*vconx_traj + gxy[i3D]*vcony_traj + gxz[i3D]*vconz_traj;
	  CCTK_REAL vcovy_traj =  gxy[i3D]*vconx_traj + gyy[i3D]*vcony_traj + gyz[i3D]*vconz_traj;
	  CCTK_REAL vcovz_traj =  gxz[i3D]*vconx_traj + gyz[i3D]*vcony_traj + gzz[i3D]*vconz_traj;

	  CCTK_REAL vsq_traj = vcony_traj*vcovx_traj + vcony_traj*vcovy_traj + vconz_traj*vcovz_traj + 1e-40;
          CCTK_REAL wlorentz_traj = sqrt(1/(1-vsq_traj));

          CCTK_REAL betaxlow = gxx[i3D]*betax[i3D] + gxy[i3D]*betay[i3D] + gxz[i3D]*betaz[i3D];
          CCTK_REAL betaylow = gxy[i3D]*betax[i3D] + gyy[i3D]*betay[i3D] + gyz[i3D]*betaz[i3D];
          CCTK_REAL betazlow = gxz[i3D]*betax[i3D] + gyz[i3D]*betay[i3D] + gzz[i3D]*betaz[i3D];
          CCTK_REAL beta2 = betax[i3D]*betaxlow + betay[i3D]*betaylow + betaz[i3D]*betazlow;
         
          CCTK_REAL rho_1, rho_2;
          if ( CCTK_EQUALS(FakeTmunu_matter, "gaussian") ) {
//FakeTmunu_matter == "gaussian") {

                CCTK_REAL centered_rad_1 = (x[i3D]-xtraj)*(x[i3D]-xtraj) +  (y[i3D]-ytraj)*(y[i3D]-ytraj) +  (z[i3D]-ztraj)*(z[i3D]-ztraj); 
                CCTK_REAL centered_rad_2 = (x[i3D]+xtraj)*(x[i3D]+xtraj) +  (y[i3D]+ytraj)*(y[i3D]+ytraj) +  (z[i3D]+ztraj)*(z[i3D]+ztraj);
 
 		rho_1 = FakeTmunu_rho * exp( - (centered_rad_1) * FakeTmunu_median ); 
 		rho_2 = FakeTmunu_rho * exp( - (centered_rad_2) * FakeTmunu_median );
          }
 
	 CCTK_REAL press_1 = initial_k*pow(rho_1, initial_Gamma); 
         CCTK_REAL eps_1 = press_1 / ((rho_1 + 1e-300) * (initial_Gamma - 1.0));

	 CCTK_REAL press_2 = initial_k*pow(rho_2, initial_Gamma); 
         CCTK_REAL eps_2 = press_2 / ((rho_2 + 1e-300) * (initial_Gamma - 1.0));

         CCTK_REAL rhoenthalpy_1 = wlorentz_traj*wlorentz_traj * (rho_1*(1.0 + eps_1) + press_1);
         CCTK_REAL rhoenthalpy_2 = wlorentz_traj*wlorentz_traj * (rho_2*(1.0 + eps_2) + press_2);


         CCTK_REAL utlow_1 = (-alp[i3D] + vconx_traj*betaxlow + vcony_traj*betaylow + vconz_traj*betazlow);
         CCTK_REAL uxlow_1 = vcovx_traj; //velxlow
         CCTK_REAL uylow_1 = vcovy_traj; //velylow
         CCTK_REAL uzlow_1 = vcovz_traj; //velzlow
         eTtt[i3D] += (rhoenthalpy_1*utlow_1*utlow_1 + press_1*(beta2 - alp[i3D]*alp[i3D]));
         eTtx[i3D] += (rhoenthalpy_1*utlow_1*uxlow_1 + press_1*betaxlow);
         eTty[i3D] += (rhoenthalpy_1*utlow_1*uylow_1 + press_1*betaylow);
         eTtz[i3D] += (rhoenthalpy_1*utlow_1*uzlow_1 + press_1*betazlow);
         eTxx[i3D] += (rhoenthalpy_1*uxlow_1*uxlow_1 + press_1*gxx[i3D]);
         eTyy[i3D] += (rhoenthalpy_1*uylow_1*uylow_1 + press_1*gyy[i3D]);
         eTzz[i3D] += (rhoenthalpy_1*uzlow_1*uzlow_1 + press_1*gzz[i3D]);
         eTxy[i3D] += (rhoenthalpy_1*uxlow_1*uylow_1 + press_1*gxy[i3D]);
         eTxz[i3D] += (rhoenthalpy_1*uxlow_1*uzlow_1 + press_1*gxz[i3D]);
         eTyz[i3D] += (rhoenthalpy_1*uylow_1*uzlow_1 + press_1*gyz[i3D]);

          
         CCTK_REAL utlow_2 = (-alp[i3D] - vconx_traj*betaxlow - vcony_traj*betaylow - vconz_traj*betazlow);
         CCTK_REAL uxlow_2 = -1.0*vcovx_traj; //velxlow
         CCTK_REAL uylow_2 = -1.0*vcovy_traj; //velylow
         CCTK_REAL uzlow_2 = -1.0*vcovz_traj; //velzlow
         eTtt[i3D] += (rhoenthalpy_2*utlow_2*utlow_2 + press_2*(beta2 - alp[i3D]*alp[i3D]));
         eTtx[i3D] += (rhoenthalpy_2*utlow_2*uxlow_2 + press_2*betaxlow);
         eTty[i3D] += (rhoenthalpy_2*utlow_2*uylow_2 + press_2*betaylow);
         eTtz[i3D] += (rhoenthalpy_2*utlow_2*uzlow_2 + press_2*betazlow);
         eTxx[i3D] += (rhoenthalpy_2*uxlow_2*uxlow_2 + press_2*gxx[i3D]);
         eTyy[i3D] += (rhoenthalpy_2*uylow_2*uylow_2 + press_2*gyy[i3D]);
         eTzz[i3D] += (rhoenthalpy_2*uzlow_2*uzlow_2 + press_2*gzz[i3D]);
         eTxy[i3D] += (rhoenthalpy_2*uxlow_2*uylow_2 + press_2*gxy[i3D]);
         eTxz[i3D] += (rhoenthalpy_2*uxlow_2*uzlow_2 + press_2*gxz[i3D]);
         eTyz[i3D] += (rhoenthalpy_2*uylow_2*uzlow_2 + press_2*gyz[i3D]);

        } 
      }
    }


}
