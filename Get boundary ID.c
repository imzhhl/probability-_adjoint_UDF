#include "udf.h"
#include "SuperUdfExtension.h"
#pragma comment(lib, "SuperUdfExtension.lib")
DEFINE_ON_DEMAND(GetOutletId)
{
	int Zone_ID;	 /* Auto get zone Id by zone name*/
	face_t f;
	cell_t c, c0;
	Thread *ft, *ct, *t0;
	Domain *domain=Get_Domain(1); /* Get the domain using Fluent utility */

#if !RP_NODE
	Zone_ID=SuperUdf_GetZoneIdByName("top"); /*get the id of zone whose name is "outlet"*/
#endif
	host_to_node_int_1(Zone_ID);
#if !RP_HOST
	/* Initialise Cells */
	/* this loops over all cells and lets the UDM = 0 */
	thread_loop_c(ct, domain)
	{
		begin_c_loop(c,ct)
		{
			C_UDMI(c,ct,0) = 0.0;
		}
		end_c_loop(c,ct)
	}
	
	/* Loop over all faces on wall */
	ft = Lookup_Thread(domain, Zone_ID);
	begin_f_loop(f,ft)
	{
		/* c0 and t0 identify the adjacent cell */
		c0 = F_C0(f, ft);
		t0 = THREAD_T0(ft);
		/* this loops over all cells adjacent to wall and lets the UDM = 2.0 */
		C_UDMI(c0, t0, 0) = 2.0;
	}
	end_f_loop(f,ft)
#endif
}

DEFINE_EXECUTE_ON_LOADING(load,libudf)
{
	SuperUdf_Initialize(AfxGetInstanceHandle());
}

