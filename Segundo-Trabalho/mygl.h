 #ifndef _MYGL_H_
#define _MYGL_H_

#include "definitions.h"
#include <cmath>
#include <iomanip>
#include <glm/glm.hpp>
#include <glm/vec4.hpp>
#include <glm/gtc/type_ptr.hpp>
//*****************************************************************************
// Defina aqui as suas funções gráficas
//*****************************************************************************

struct Color
{
	unsigned char r, g, b, a;

	Color(unsigned char _r, unsigned char _g, unsigned char _b, unsigned char _a){
		r = _r;
		g = _g;
		b = _b;
		a = _a;
	}
};

class Vertex
{
public:

	unsigned int x, y, z;
	unsigned char r, g, b, a;


	Vertex(unsigned int x, unsigned int y, unsigned int z, unsigned char r, unsigned char g, unsigned char b, unsigned char a){
		this->x = x;
		this->y = y;
		this->z = z;
		this->r = r;
		this->g = g;
		this->b = b;
		this->a = a;
	}

	Vertex(unsigned int x, unsigned int y, unsigned char r, unsigned char g, unsigned char b, unsigned char a){
		this->x = x;
		this->y = y;
		this->z = 0;
		this->r = r;
		this->g = g;
		this->b = b;
		this->a = a;
	}

	Vertex(unsigned int x, unsigned int y, unsigned int z, const Color& color){
		this->x = x;
		this->y = y;
		this->z = z;
		r = color.r;
		g = color.g;
		b = color.b;
		a = color.a;
	}

	Vertex(glm::vec4 v, Color& color){
		x = v.x;
		y = v.y;
		z = v.z;
		r = color.r;
		g = color.g;
		b = color.b;
		a = color.a;
	}
};

void PutPixel(const Vertex vertex){
int offset = 4*vertex.x + 4*IMAGE_WIDTH*vertex.y;
	double pctAlpha = vertex.a/255.0;
	if(vertex.x <=511 && vertex.y<=511 && vertex.x>=0 && vertex.y>=0){
		FBptr[offset] = vertex.r*(pctAlpha);
		FBptr[offset + 1] = vertex.g*(pctAlpha);
		FBptr[offset + 2] = vertex.b*(pctAlpha);
		FBptr[offset + 3] = vertex.a;
	}

}

void DrawLine(const Vertex vertex1, const Vertex vertex2){
	
	if(vertex1.x == vertex2.x && vertex1.y == vertex2.y){
		PutPixel(vertex1);
		return;
	}

	int dx = vertex2.x - vertex1.x;
	int dy = vertex2.y - vertex1.y;
	double m = (double)dy/dx;
	double dist = sqrt(pow(dx,2)+ pow(dy,2));
	double difR = ((int)vertex2.r - (int)vertex1.r)/dist;
	double difG = ((int)vertex2.g - (int)vertex1.g)/dist;
	double difB = ((int)vertex2.b - (int)vertex1.b)/dist;
	double difA = ((int)vertex2.a - (int)vertex1.a)/dist;

	int i=0;

//--------------------------------------------------	
	if(!dx){
		if(dy > 0){
			for(int i =0; i <= dy; i++){
				PutPixel(Vertex(vertex1.x, vertex2.y-i, vertex2.r - difR*i , vertex2.g - difG*i , vertex2.b - difB*i,  vertex2.a - difA*i));
			}
		}else{
			for(int i =0; i >= dy; i--){
				PutPixel(Vertex(vertex1.x, vertex1.y+i, vertex1.r - difR*i , vertex1.g - difG*i , vertex1.b - difB*i,  vertex1.a - difA*i));
			}
		}
		return;
	}
//--------------------------------------------------	
	else if(m>0 && m<1){
		int dif = (dy>0)? 1:-1;
		const Vertex *v = (dy>0)? &vertex1:&vertex1;
		//std::cout << dif << std::endl;
		int d = (2 * dy - dx) * dif;
		int incr_e = 2 * dy * dif;
		int incr_ne = 2 * (dy - dx) * dif;
		int x ,y;
		x = vertex1.x;
		y = vertex1.y;
		PutPixel(Vertex(x, y, vertex1.r, vertex1.g, vertex1.b, 1));
		while ((dy>0)? (x < vertex2.x):(x > vertex2.x)) {			
			if (d <= 0) {
				d += incr_e;
				x += dif;
			} else {
				d += incr_ne;
				x += dif;
				y += dif;				
			}
			i++;
			PutPixel(Vertex(x, y, v->r + difR*i , v->g + difG*i , v->b + difB*i, v->a + difA*i));
		}
	}
//--------------------------------------------------	
	else if(m<0 && m>-1){
		int dif = (dx>0)? 1:-1;
		const Vertex *v = (dy>0)? &vertex1:&vertex1;
		//std::cout << dif << std::endl;
		int d = (2 * dy + dx) * dif;
		int incr_e = -2 * dy * dif;
		int incr_ne = -2 * (dy + dx) * dif;
		//std::cout << "d: " << d << " incr_e: " << incr_e << " incr_ne: "  << incr_ne<< std::endl;
		int x ,y;
		x = vertex1.x;
		y = vertex1.y;
		PutPixel(Vertex(x, y, vertex1.r, vertex1.g, vertex1.b, 1));
		while ((dy>0)? (x > vertex2.x):(x < vertex2.x)) {			
			if (d <= 0) {

				d += incr_e;
				x += dif;
			} else {
				d += incr_ne;
				x += dif;
				y -= dif;				
			}
			i++;
			PutPixel(Vertex(x, y, v->r + difR*i , v->g + difG*i , v->b + difB*i, v->a + difA*i));
		}
	}
//--------------------------------------------------	
	else if(m>1){
		int dif = (dy>0)? 1:-1;
		const Vertex *v = (dy>0)? &vertex1:&vertex1;
		int d = (2 * dx - dy) * dif;
		int incr_e = 2 * dx * dif;
		int incr_ne = 2 * (dx - dy) * dif;
		int x = vertex1.x;
		int y = vertex1.y;
		PutPixel(Vertex(x, y, vertex1.r, vertex1.g, vertex1.b, vertex1.a));
		while ((dy>0)? (y < vertex2.y):(y > vertex2.y)) {
			if (d <= 0) {
				d += incr_e;
				y+=dif;
			} else {
				d += incr_ne;
				x+=dif;
				y+=dif;
			}
			++i;
			PutPixel(Vertex(x, y, v->r + difR*i , v->g + difG*i , v->b + difB*i, v->a + difA*i));
		}
	}	
//--------------------------------------------------
	else if(m<-1){
		int dif = (dy<0)? 1:-1;
		const Vertex *v = (dy<0)? &vertex1:&vertex1;
		int d = (2 * dx - dy) * dif;
		int incr_e = 2 * dx * dif;
		int incr_ne = 2 * (dx + dy) * dif;
		int x = vertex1.x;
		int y = vertex1.y;
		PutPixel(Vertex(x, y, vertex1.r, vertex1.g, vertex1.b, vertex1.a));
		while ((dy>0)? (y < vertex2.y):(y > vertex2.y)) {
			if (d <= 0) {
				d += incr_e;
				y+=-dif;
			} else {
				d += incr_ne;
				x+=dif;
				y+=-dif;
			}
			++i;
			PutPixel(Vertex(x, y, v->r + difR*i , v->g + difG*i , v->b + difB*i, v->a + difA*i));
		}
	}
//--------------------------------------------------
	else if(!m){
		if(dx > 0){
			for(int i =0; i <= dx; i++){
				PutPixel(Vertex(vertex1.x+i, vertex1.y, vertex1.r + difR*i , vertex1.g + difG*i , vertex1.b + difB*i,  vertex1.a + difA*i));
			}
		}else{
			for(int i =0; i >= dx; i--){
				PutPixel(Vertex(vertex1.x+i, vertex1.y, vertex1.r - difR*i , vertex1.g - difG*i , vertex1.b - difB*i,  vertex1.a - difA*i));
			}
		}
	}
//--------------------------------------------------
	else if(m == 1){
		if(dy > 0){
			for(int i =0; i <= dy; i++){
				PutPixel(Vertex(vertex1.x+i, vertex1.y+i, vertex1.r + difR*i , vertex1.g + difG*i , vertex1.b + difB*i,  vertex1.a + difA*i));
			}
		}else{
			for(int i =0; i >= dy; i--){
				PutPixel(Vertex(vertex1.x+i, vertex1.y+i, vertex1.r - difR*	i , vertex1.g - difG*i , vertex1.b - difB*i,  vertex1.a - difA*i));
			}
		}
	}
//--------------------------------------------------
	else{
		if(dy > 0){
			for(int i =0; i <= dy; i++){
				PutPixel(Vertex(vertex1.x-i, vertex1.y+i, vertex1.r + difR*i , vertex1.g + difG*i , vertex1.b + difB*i,  vertex1.a + difA*i));
			}
		}else{
			for(int i =0; i >= dy; i--){
				PutPixel(Vertex(vertex1.x-i, vertex1.y+i, vertex2.r - difR*i , vertex2.g - difG*i , vertex2.b - difB*i,  vertex1.a + difA*i)); 	
			}
		}
	}
}

void DrawTriangle(const Vertex vertex1, const Vertex vertex2, const Vertex vertex3)
{
		DrawLine(vertex1, vertex2);
		DrawLine(vertex2, vertex3);
		DrawLine(vertex1, vertex3);	
}


glm::mat4 getModel(const glm::vec3& pos = glm::vec3(), const glm::vec3& rot = glm::vec3(), const glm::vec3& scale = glm::vec3(1.0f, 1.0f, 1.0f))
{
    glm::mat4 posMatrix(1,0,0,pos.x,
    					0,1,0,pos.y,
    					0,0,1,pos.z,
    					0,0,0,1);

    glm::mat4 rotXMatrix(1,0		 ,0			 ,0,
    					 0,cos(rot.x),-sin(rot.x),0,
    					 0,sin(rot.x), cos(rot.x),0,
    					 0,0          ,0         ,1);
    glm::mat4 rotYMatrix(cos(rot.y) ,0,sin(rot.y),0,
    					 0          ,1,0         ,0,
    					 -sin(rot.y),0,cos(rot.y),0,
    					 0          ,0,0         ,1);
    glm::mat4 rotZMatrix(cos(rot.z),-sin(rot.z),0        ,0,
    					 sin(rot.z),cos(rot.z) ,0        ,0,
    					 0         ,0         ,1         ,0,
    					 0         ,0         ,0         ,1);
    glm::mat4 scaleMatrix(scale.x,0,0,0,
    					  0,scale.y,0,0,
    					  0,0,scale.z,0,
    					  0,0,0,1);

    glm::mat4 rotMatrix = rotZMatrix * rotYMatrix * rotXMatrix;

    return  posMatrix * rotMatrix  * scaleMatrix ;
}

glm::mat4 getView(const glm::vec3& pos, const glm::vec3& lookat, const glm::vec3& up)
{

	glm::vec3 z_camera = (pos - lookat)/glm::length(pos - lookat);
	glm::vec3 x_camera = glm::cross(up, z_camera) / glm::length(glm::cross(up, z_camera));
	glm::vec3 y_camera = glm::cross(z_camera, x_camera) / glm::length(glm::	cross(z_camera, x_camera));

	glm::mat4 Bt(x_camera.x, x_camera.y, x_camera.z, 0, 
			     y_camera.x, y_camera.y, y_camera.z, 0, 
			     z_camera.x, z_camera.y, z_camera.z, 0, 
			     0		   , 0		   , 0		   , 1);
	glm::mat4 T(1, 0, 0, -pos.x, 
		        0, 1, 0, -pos.y, 
		        0, 0, 1, -pos.z, 
		        0, 0, 0,      1);
	return Bt * T;	
}

glm::mat4 getProjection(double d)
{
	glm::mat4 projection(1.0,0.0,0.0,0.0,
                		 0.0,1.0,0.0,0.0,
                		 0.0,0.0,1.0,d,
                		 0.0,0.0,-1/d,1);
	return projection;
}
/*
glm::mat4 getProjection(float fov, float aspect, float near, float far)
{
	float aux = (-near+far)/near;
	float aux2 = -(1.0/near);
	glm::mat4 projection(1.0,0.0,0.0,0.0,
                		 0.0,1.0,0.0,0.0,
                		 0.0,0.0,aux,far,
                		 0.0,0.0,aux2,1.0);
	return projection;
}*/

glm::mat4 getViewport(){
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
	return s2*T*s1;
}





#endif // _MYGL_H_

