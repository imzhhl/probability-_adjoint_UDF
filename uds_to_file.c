/******************************************************************************
* Created on : 2022-09-29 17:12
* author : ZHHL_WHU
* 功能：将对偶uds场存入文件中，然后转存到下一个case中，进行逆时间迭代
******************************************************************************/

#include "udf.h"
#define FLUID_ID 5  /*需要转存uds的计算域的ID*/

DEFINE_ON_DEMAND(uds_to_file)
{
	#if  !RP_HOST
		Domain* domain = Get_Domain(1); /*获取计算域*/
		Thread* t; /*获取指针*/
		cell_t c;   /*标定单元*/
	#else
		int i;  /*用来循环计数*/
	#endif

	#if !RP_NODE
		FILE *fp = NULL;   /*定义一个文件指针*/
		char filename[] = "uds.txt";   /*文件名*/
	#endif

	#if PARALLEL
		int size;   /*传递变量的数据,实际上计算节点上的单元数*/
		real* uds_1;    /*uds_1的值*/
		real* uds_2;    /*uds_1的值*/
		real* uds_3;    /*uds_1的值*/
		int pe; /*用来保存主机和计算节点*/
	#endif

	#if !RP_HOST
		t = Lookup_Thread(domain, FLUID_ID);
	#endif

	#if !RP_NODE
		if ((fp = fopen(filename, "w")) == NULL)   /*判断文件是否打开成功*/
			Message("\n Warning: Unable to open %s for writing\n", filename);   /*文件打开失败*/
		else
			Message("\nWriting uds to %s...", filename);   /*打开文件成功*/
	#endif

	#if !PARALLEL   /*在串行情况下，输出文件*/
		begin_c_loop(c, t)
			fprintf(fp, "%g, %g, %g\n", C_UDSI(c, t, 1), C_UDSI(c, t, 2), C_UDSI(c, t, 3));
		end_c_loop(c,t)
	#endif

	#if RP_NODE
		size=THREAD_N_ELEMENTS_INT(t);  /*得到分块计算域内单元的数目*/
		uds_1 = (real *)malloc(size * sizeof(real)); /*为uds_1分配空间*/
		uds_2 = (real *)malloc(size * sizeof(real)); /*为uds_2分配空间*/		
		uds_3 = (real *)malloc(size * sizeof(real)); /*为uds_3分配空间*/		
		begin_c_loop_int(c, t)
		{
			uds_1[c] = C_UDSI(c, t, 1);   /*得到uds_1的值*/
			uds_2[c] = C_UDSI(c, t, 2);   /*得到uds_2的值*/
			uds_3[c] = C_UDSI(c, t, 3);   /*得到uds_3的值*/
		}
		end_c_loop_int(c, t)

		/***************************************
		*判断是否是0计算节点
		*当前节点为0计算节点，则pe为node_host,
		*通过0计算节点向主机节点发送数据
		*当前节点为非0计算节点，则pe为node_zero
		*其他计算节点向0计算节点发送数据
		****************************************/
		pe = (I_AM_NODE_ZERO_P) ? node_host: node_zero;
		/********************************************
		*当前节点为0节点，则当前节点向主机发送数据
		*当前节点为非0节点，则当前节点向0节点发送数据
		*********************************************/

		PRF_CSEND_INT(pe, &size, 1, myid);
		PRF_CSEND_REAL(pe, uds_1, size, myid);
		PRF_CSEND_REAL(pe, uds_2, size, myid);
		PRF_CSEND_REAL(pe, uds_3, size, myid);

		/*发送数据以后，释放内存*/
		free(uds_1);
		free(uds_2);
		free(uds_3);

		if (I_AM_NODE_ZERO_P)   /*如果是0计算节点*/
		{
			compute_node_loop_not_zero (pe)  /*循环非0节点*/
			{
				PRF_CRECV_INT(pe, &size, 1, pe);   /*0计算节点从其他计算节点接收各分块节点的单元数*/
				uds_1 = (real *)malloc(size * sizeof(real)); /*为uds_1分配空间*/
				uds_2 = (real *)malloc(size * sizeof(real)); /*为uds_1分配空间*/
				uds_3 = (real *)malloc(size * sizeof(real)); /*为uds_1分配空间*/

				/*0计算节点从其他计算节点接收各结算节点的单元对应的变量*/
				PRF_CRECV_REAL(pe, uds_1, size, pe);
				PRF_CRECV_REAL(pe, uds_2, size, pe);
				PRF_CRECV_REAL(pe, uds_3, size, pe);
				
				/*0计算节点将从各分块节点单元数发送给主机节点*/
				PRF_CSEND_INT(node_host, &size, 1, myid);   

				/*0计算节点将从其他计算节点接收各结算节点的单元对应的变量发送给主机节点*/
				PRF_CSEND_REAL(node_host, uds_1, size, myid);  
				PRF_CSEND_REAL(node_host, uds_2, size, myid); 
				PRF_CSEND_REAL(node_host, uds_3, size, myid); 			
				
				/*释放内存*/
				free((char *)uds_1);
				free((char *)uds_2);
				free((char *)uds_3);
			}
		}
	#endif

	#if RP_HOST
		compute_node_loop(pe)   /*仅仅作为一个循环计数器的角色*/
		{
			/*接收各节点的数据*/
			PRF_CRECV_INT(node_zero, &size, 1, node_zero);
			uds_1 = (real *)malloc(size * sizeof(real));
			uds_2 = (real *)malloc(size * sizeof(real));
			uds_3 = (real *)malloc(size * sizeof(real));

			PRF_CRECV_REAL(node_zero, uds_1, size, node_zero);  
			PRF_CRECV_REAL(node_zero, uds_2, size, node_zero);  
			PRF_CRECV_REAL(node_zero, uds_3, size, node_zero);  

			/*将数据写入文件当中*/

			for (i=0; i<size; i++)
			{
				fprintf(fp, "%g,%g,%g\n", uds_1[i], uds_2[i], uds_3[i]);
			}
			/*释放内存*/
			free(uds_1);
			free(uds_2);
			free(uds_3);
	
		}
	#endif

	#if !RP_NODE
		fclose(fp); /*关闭文件*/
		Message("Done\n");  /*输出关闭文件完成*/
	#endif
}
