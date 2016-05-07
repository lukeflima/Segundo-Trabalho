#include "main.h"
#include "mygl.h"
#include <cmath>
#include <vector>
#include "objLoader.h"
#define PI 3.14159265
/*
void DrawSquare(int x,int y, int size, int i){
	Vertex v1(x*sin(i*PI/180.0),y,255,255,255,255);
	Vertex v2((x+size)*sin(i*PI/180.0), y, 255, 255, 255, 255);
	Vertex v3((x+size)*sin(i*PI/180.0), y+size, 255,255,255,255);
	Vertex v4((x)*sin(i*PI/180.0), y+size, 255, 255, 255,255);
	DrawLine(v1,v2);
	DrawLine(v2,v3);
	DrawLine(v3,v4);
	DrawLine(v4,v1);
}
*/
void Clear(){
	for(int i =0; i < 4*IMAGE_WIDTH*IMAGE_HEIGHT; i++)
	{
		FBptr[i] = 0;
		FBptr[i+1] = 0;
		FBptr[i+2] = 0;
		FBptr[i+3] = 0;
	}
}
objLoader *obj;
//-----------------------------------------------------------------------------
void MyGlDraw(void)
{
	//*************************************************************************
	// Chame aqui as funções do mygl.h
	//*************************************************************************	
	Clear();
	static int a = 0;
	//double sina = cos(0.01*a);
	
	obj->load("monkey_head2.obj");

	Color color(255,255,255,255);

	std::vector<glm::vec4> list;
	glm::vec4 x(1,0,0,1);
	glm::vec4 y(0,1,0,1);
	glm::vec4 z(0,0,1,1);

	glm::mat4 modelviewprojection = getProjection(5)* 
				  (getView(glm::vec3(0,0,5),  glm::vec3(0,0,0), glm::vec3(0,1,0))*
				  getModel(glm::vec3(0,0,0),  glm::vec3(0,0.01*a,0),  glm::vec3(0.6f,0.6f,0.6f)));
	/*	const float *pSource = (const float*)glm::value_ptr(modelviewprojection);
	std::cout << pSource[0] << ',' << pSource[1] << ',' << pSource[2] << ',' << pSource[3] << std::endl
			  << pSource[4] << ',' << pSource[5] << ',' << pSource[6] << ',' << pSource[7] << std::endl
			  << pSource[8] << ',' << pSource[9] << ',' << pSource[10] << ',' << pSource[11] << std::endl
			  << pSource[12] << ',' << pSource[13] << ',' << pSource[14] << ',' << pSource[15] << std::endl  << std::endl;
	*/

	glm::mat4 viewport = getViewport();
	
		glm::mat4 s1(1,0,0,0,
			   		 0,-1,0,0,
			   		 0,0,1,0,
			    	 0,0,0,1);
	glm::mat4 T(1,0,0,1,
			    0,1,0,1,
			    0,0,1,0,
			    0,0,0,1);
	glm::mat4 s2(511/2.0,     0,0,0,
				     0, 511/2.0,0,0,
				     0,     0,1,0,
				     0,     0,0,1);

	for (int i = 0; i < obj->vertexCount; ++i)
	{	
		float x = obj->vertexList[i]->e[0],y = obj->vertexList[i]->e[1],z = obj->vertexList[i]->e[2], w = 1.0f;
		list.push_back(glm::vec4(x,y,z,w));

		list[i]= modelviewprojection*list[i];
		
		list[i].x /= list[i][3];
		list[i].y /= list[i][3];
		list[i].z /= list[i][3];
		list[i][3] /= list[i][3];

		
		list[i]= s1*list[i];		
 		list[i]= list[i]*T;	
		list[i]= s2*list[i];		
		
	}

			
	
	for(unsigned int i=0; i<obj->faceCount; i++){
		obj_face *obj2 = obj->faceList[i];
		if(list[obj2->vertex_index[0]].x>511 || list[obj2->vertex_index[1]].x>511||list[obj2->vertex_index[2]].x> 511||
			list[obj2->vertex_index[0]].x<0 || list[obj2->vertex_index[1]].x<0||list[obj2->vertex_index[2]].x< 0||
			list[obj2->vertex_index[0]].y>511 || list[obj2->vertex_index[1]].y>511||list[obj2->vertex_index[2]].y> 511||
			list[obj2->vertex_index[0]].y<0 || list[obj2->vertex_index[1]].y<0||list[obj2->vertex_index[2]].y< 0);
		else DrawTriangle(Vertex(list[obj2->vertex_index[0]], color),
					 Vertex(list[obj2->vertex_index[1]], color),
					 Vertex(list[obj2->vertex_index[2]], color));
		
	}	
	modelviewprojection = getProjection(5)* 
				  (getView(glm::vec3(0,0,5),  glm::vec3(0,0,0), glm::vec3(0,1,0))*
				  getModel(glm::vec3(0,0,0),  glm::vec3(0,0.01*a,0),  glm::vec3(1.0f,1.0f,1.0f)));
	
	x= modelviewprojection*x;
		
		x.x /= x[3];
		x.y /= x[3];
		x.z /= x[3];
		x[3] /= x[3];

		
		x= s1*x;		
 		x= x*T;	
		x= s2*x;		
	y= modelviewprojection*y;
		
		y.x /= y[3];
		y.y /= y[3];
		y.z /= y[3];
		y[3] /= y[3];

		
		y= s1*y;		
 		y= y*T;	
		y= s2*y;
		z= modelviewprojection*z;
		
		z.x /= z[3];
		z.y /= z[3];
		z.z /= z[3];
		z[3] /= z[3];

		
		z= s1*z;		
 		z= z*T;	
		z= s2*z;

	Color cx(255,0,0,255);
	Color cy(0,255,0,255);
	Color cz(0,0,255,255);
	DrawLine(Vertex(glm::vec4(255,255,0,1),cx), Vertex(x, cx));
	DrawLine(Vertex(glm::vec4(255,255,0,1),cy), Vertex(y, cy));
	DrawLine(Vertex(glm::vec4(255,255,0,1),cz), Vertex(z, cz));
	a++;
	/*
	Clear();
	static int i;
	
	DrawSquare(200,200,100,i);
	DrawSquare(125,125,100,i);
	DrawSquare(275,125,100,i);
	DrawSquare(200,50,100,i);

	DrawLine(Vertex(125*sin(i*PI/180.0),125,255,255,255,255), Vertex(200*sin(i*PI/180.0),50,255,255,255,255));
	DrawLine(Vertex(125*sin(i*PI/180.0),125+100,255,255,255,255), Vertex(200*sin(i*PI/180.0),300,255,255,255,255));
	DrawLine(Vertex(300*sin(i*PI/180.0),50,255,255,255,255), Vertex((125+250)*sin(i*PI/180.0),125,255,255,255,255));
	DrawLine(Vertex((125+250)*sin(i*PI/180.0),125+100,255,255,255,255), Vertex(300*sin(i*PI/180.0),300,255,255,255,255));

	DrawLine(Vertex((125+100)*sin(i*PI/180.0),125,255,255,255,255), Vertex((200+100)*sin(i*PI/180.0),50,255,255,255,255));
	DrawLine(Vertex((125+100)*sin(i*PI/180.0),125+100,255,255,255,255), Vertex((200+100)*sin(i*PI/180.0),300,255,255,255,255));
	DrawLine(Vertex(300*sin(i*PI/180.0),50+100,255,255,255,255), Vertex((125+250)*sin(i*PI/180.0),125+100,255,255,255,255));
	DrawLine(Vertex((125+250)*sin(i*PI/180.0),125+100-100,255,255,255,255), Vertex(300*sin(i*PI/180.0),300-100,255,255,255,255));

	DrawLine(Vertex(125*sin(i*PI/180.0),125+100,255,255,255,255), Vertex(200*sin(i*PI/180.0),50+100,255,255,255,255));
	DrawLine(Vertex(125*sin(i*PI/180.0),125,255,255,255,255), Vertex(200*sin(i*PI/180.0),300-100,255,255,255,255));
	DrawLine(Vertex((300-100)*sin(i*PI/180.0),50,255,255,255,255), Vertex((125+250-100)*sin(i*PI/180.0),125,255,255,255,255));
	DrawLine(Vertex((125+250-100)*sin(i*PI/180.0),125+100,255,255,255,255), Vertex((300-100)*sin(i*PI/180.0),300,255,255,255,255));

	DrawLine(Vertex(275*sin(i*PI/180.0),125,255,255,255,255), Vertex(200*sin(i*PI/180.0),200,255,255,255,255));
	DrawLine(Vertex(300*sin(i*PI/180.0),150,255,255,255,255), Vertex(225*sin(i*PI/180.0),225,255,255,255,255));
	DrawLine(Vertex(225*sin(i*PI/180.0),125,255,255,255,255), Vertex(300*sin(i*PI/180.0),200,255,255,255,255));
	DrawLine(Vertex(200*sin(i*PI/180.0),150,255,255,255,255), Vertex(275*sin(i*PI/180.0),225,255,255,255,255));
	
	i+= 1;

	*/

}

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
	obj = new objLoader();
	// Inicializações.
	InitOpenGL(&argc, argv);
	InitCallBacks();
	InitDataStructures();

	// Ajusta a função que chama as funções do mygl.h
	DrawFunc = MyGlDraw;	

	// Framebuffer scan loop.
	glutMainLoop();

	return 0;
}

