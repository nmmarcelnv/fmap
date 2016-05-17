#include "zdread.h"

int ZDread(char *zdockfile, ZDCOM *zdcom, double angl[][3], double tran[][3], double zdscore[], bool load)
{
	FILE *fpZD;
	int nLine,nPose;
	char line[MAXLENLINE];
	double a1,a2,a3,t1,t2,t3,score;
	if ((fpZD=fopen(zdockfile, "r"))==NULL) 
	{
		printf("Cannot open file:%s\n",zdockfile);
		exit(1);
	}
	nLine=0;
	nPose=0;
	while (fgets(line, MAXLENLINE, fpZD) != NULL)
	{
		switch (nLine)
		{
			case 0:
			{
				sscanf(line,"%d%lf",&(zdcom->N),&(zdcom->spacing));
				//printf("%-.0f\t%-.1f\n",zdcom->N,zdcom->spacing);
				break;
			}
			case 1:
			{
			
				sscanf(line,"%lf%lf%lf",&(zdcom->rand[0]),&(zdcom->rand[1]),&(zdcom->rand[2]));
				//printf("%lf\t%lf\t%lf\n",zdcom->rand[0],zdcom->rand[1],zdcom->rand[2]);
				break;
			}
			case 2:
			{
				sscanf(line,"%s%lf%lf%lf",(char *)&(zdcom->rec),&(zdcom->r[0]),&(zdcom->r[1]),&(zdcom->r[2]));
				//printf("%s\t%-.3f\t%-.3f\t%-.3f\n",zdcom->rec,zdcom->r[0],zdcom->r[1],zdcom->r[2]);
				break;
			}
			case 3:
			{
				sscanf(line,"%s%lf%lf%lf",(char *)&(zdcom->lig),&(zdcom->l[0]),&(zdcom->l[1]),&(zdcom->l[2]));
				//printf("%s\t%-.3f\t%-.3f\t%-.3f\n",zdcom->lig,zdcom->l[0],zdcom->l[1],zdcom->l[2]);
				break;
			}
			default:
			{
				if (load){
					sscanf(line,"%lf%lf%lf%lf%lf%lf%lf",&a1,&a2,&a3,&t1,&t2,&t3,&score);
					//printf("%lf\t%lf\t%lf\t%-.0f\t%-.0f\t%-.0f\t%6.2f\n",a1,a2,a3,t1,t2,t3,score);
					angl[nPose][0]=a1;
					angl[nPose][1]=a2;
					angl[nPose][2]=a3;
					tran[nPose][0]=t1;
					tran[nPose][1]=t2;
					tran[nPose][2]=t3;
					zdscore[nPose]=score;
				}
				nPose++;	
				break;	
			}
		}	
		nLine++;
	}
	fclose(fpZD);
	return nPose;
}
