
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

// define the base addresses
#define VGA_BUFFER 0xff203020
#define FPGA_CHAR_BASE 0x09000000
#define TIMERSEC 10000000
#define PS2_BASE 0xFF200100
#define LEDS_BASE 0xFF200000
#define HEX3_HEX0_BASE 0xFF200020
#define HEX5_HEX4_BASE 0xFF200030
#define action1 (1 << 0)
#define action2 ( 1 << 1)

// define constants
#define XRES 320
#define YRES 240
#define G 1E-5
#define SUN_MASS 1E5
#define EARTH_MASS 1
#define EARTH_VELOCITY 0.15
#define SUN_TO_EARTH 100
//0.000000000066743
#define ERRORCODE (-10000)
#define PI 3.14159
#define DOWNSCALE 1
#define XRES 320
#define YRES 240

#define BOXPIXELS 15
#define TIME_STEP 40

int NUM_PLANETS = 0;
double simulation_time = 0.0;


volatile int * pixel_ctrl_ptr = (int *)0xFF203020;
volatile int *pixel_status = (int*)0xFF20302C;
int pixel_buffer_start;

//planet pass in struct
typedef struct{
	double mass;
	double x_velocity;
	double y_velocity;
	double distanceFromSun;
    	double x_pos;
    	double y_pos;
	int planetColor;
} planetPassIn;

planetPassIn planets[5];
int planetIndex = -1;

// MOUSE STRUCTS
typedef struct PS2_Mouse{
    uint8_t status;
    int8_t x_mov;
    int8_t y_mov;

} PS2_Mouse;

typedef struct PS2{
    uint32_t data;
    uint32_t control;
} PS2;

volatile PS2* ps2 = (PS2*)PS2_BASE;

struct Point{
    uint16_t x, y;
	uint16_t prevx, prevy;
	bool toggleHoldPlanet;
};
struct Point clickPoint;

// function prototypes
void HEX_PS2(char b1, char b2, char b3);
void drawPNG(const uint16_t arr[240][320]);
void plot_pixel(int, int, uint16_t);
void clear_screen();
void video_text_delay(int x, int y, char * text_ptr);
void video_text_instant(int x, int y, char * text_ptr);
void draw_point(struct Point p, short int colour);
void wait_for_vsync();
void waitasec(int pow_fraction);
void PS2_send(uint8_t);
uint8_t PS2_get();
void swapBuffer();
void drawPNG(const uint16_t arr[YRES][XRES]);
void redraw_background_for_mouse(struct Point p, const uint16_t arr[YRES][XRES]);

struct timer_t {
       volatile unsigned int status;
       volatile unsigned int control;
       volatile unsigned int periodlo;
       volatile unsigned int periodhi;
       volatile unsigned int snaplo;
       volatile unsigned int snaphi;
};


typedef struct Kepler_Variables {

    //initial given values and orbit position varaibles
    double mass;
    double position[2];
    double velocity[2];
    double phi;
    double r_init;

    double radius;
    double angle;
    double previous_position[2];
    bool active;

    //orbital description variables
    double eccentricity;
    double Energy;
    double angular_moment;
    double semi_major_axis;
    int planetColor;

} Kepler_Variables;


double getEnergy(double m, double v, double M, double r);
double getAngularMomentum(double v, double r_init, double, double phi);
double getEccentricity(double L, double E, double m, double M);


double radius_WRT_angle(double angle, double L, double m, double M, double e);
void update_position_given_polar(double* position, double angle, double r);
void full_Kepler(struct Kepler_Variables* planet, double* sun_position, double sun_mass);
//updates orbit position based on the Kepler vars, creating a function of time: r(t)
void orbit(struct Kepler_Variables* planet, double sun_mass);
void init_planet(struct Kepler_Variables* planet, int planetColor, double mass, double* position, double* velocity, double* sun_pos);
void update_velocity(Kepler_Variables* planet, double angle, double M);


double rad_to_angle(double input);
double angle_to_rad(double input);
double dotproduct(double* v1, double* v2, int size);
void difference_vector(double* diff, double* v1, double* v2, int size);
double angle_between_2vectors(double* v1, double* v2);


void crossProduct(double v_A[], double v_B[], double cross_P[]) {
    cross_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
    cross_P[1] = v_A[2] * v_B[0] - v_A[0] * v_B[2];
    cross_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}


double rad_to_angle(double input) {
    return input * 180 / PI;
}

double angle_to_rad(double input) {
    return input * PI / 180;
}

double getEnergy(double m, double v, double M, double r) {

    return 0.5 * m * v * v - G * M * m / r;

}

double getAngularMomentum(double v, double r_init, double m, double phi) {

    return m * v * r_init * sin(phi);

}

double getEccentricity(double L, double E, double m, double M) {

    double eccentricity = (2 * L * L * E);
    eccentricity = eccentricity / (G * G * m * m * m * M * M);
    eccentricity = eccentricity + 1;
    eccentricity = sqrt(eccentricity);
    return eccentricity;

}

double getSemiMajorAxis(double E, double m, double M) {
    double a = -1 * G * M * m / 2 / E;
    return a;
}


double radius_WRT_angle(double angle, double L, double m, double M, double e) {

    double radius = (L * L / G / m / m / M / (1 + e * cos(angle)));
    return radius;

    

}


//changes the values of position passed by pointer
void update_position_given_polar(double* position, double angle, double r) {


    // angle to radians
    position[0] = r * cos(angle);
    position[1] = r * sin(angle);

}

double dotproduct(double* v1, double* v2, int size) {

    double dotproduct = 0;

    for (int i = 0; i < size; i++) {
        dotproduct += v1[i] * v2[i];
    }

    return dotproduct;
}

//+
double angle_between_2D_vectors(double* v1, double* v2) {


    double dotprod = dotproduct(v1, v2, 2);
    if (dotprod == ERRORCODE) {
        return __LINE__;
    }
    double mag1 = 0.0;
    double mag2 = 0.0;

    mag1 = sqrt((v1[0]) * (v1[0]) + (v1[1]) * (v1[1]));
    mag2 = sqrt((v2[0]) * (v2[0]) + (v2[1]) * (v2[1]));

    double angle = acos(dotprod / (mag1 * mag2));

    return angle;

}

//parameters: difference vector, v1, v2. Performs difference = v2 - v1.
void difference_vector(double* diff, double* v1, double* v2, int size) {

    for (int i = 0; i < size; i++) {
        diff[i] = v2[i] - v1[i];
    }

}

//given the two masses, initial position of both planets and velocity of small planet, the full orbit can be calculated.
void full_Kepler(struct Kepler_Variables* planet, double* sun_position, double sun_mass) {

    double v_mag = sqrt(planet->velocity[0] * planet->velocity[0] + planet->velocity[1] * planet->velocity[1]);

    planet->r_init = sqrt((planet->position[0] - sun_position[0]) * (planet->position[0] - sun_position[0]) + (planet->position[1] - sun_position[1]) * (planet->position[1] - sun_position[1]));
    double position_vector[2];
    difference_vector(position_vector, sun_position, planet->position, 2);
    planet->phi = angle_between_2D_vectors(planet->velocity, position_vector);


    planet->Energy = getEnergy(planet->mass, v_mag, sun_mass, planet->r_init);
    planet->angular_moment = getAngularMomentum(v_mag, planet->r_init, planet->mass, planet->phi);
    planet->eccentricity = getEccentricity(planet->angular_moment, planet->Energy, planet->mass, sun_mass);
    planet->semi_major_axis = getSemiMajorAxis(planet->Energy, planet->mass, sun_mass);

}


void update_velocity(Kepler_Variables* planet, double angle, double M) {

    //update magnitude of speed based on conservation of energy, then update the direction of it. 
    //v = sqrt (GM / r)
    double speed;

    double dir_x = -1 * sin(angle) / sqrt(1 + planet->eccentricity * planet->eccentricity + 2 * planet->eccentricity * cos(angle));
    double dir_y = (planet->eccentricity + cos(angle)) / sqrt(1 + planet->eccentricity * planet->eccentricity + 2 * planet->eccentricity * cos(angle));

    speed = sqrt(G * M * (2 / planet->radius - 1 / planet->semi_major_axis));

    planet->velocity[0] = speed * dir_x;
    planet->velocity[1] = speed * dir_y;

}

//updates orbit position based on the Kepler vars, as a function of time: r(time)
void orbit(struct Kepler_Variables* planet, double sun_mass) {

    //ouble v_mag = sqrt(planet->velocity[0] * planet->velocity[0] + planet->velocity[1] * planet->velocity[1]);
    double d_mag = sqrt(planet->position[0] * planet->position[0] + planet->position[1] * planet->position[1]);


    planet->radius = radius_WRT_angle(planet->angle, planet->angular_moment, planet->mass, sun_mass, planet->eccentricity);

    double position_vec[3] = { 0, 0, 0 };
    position_vec[0] = planet->position[0] / d_mag;
    position_vec[1] = planet->position[1] / d_mag; //unit vectorize the position
    double tangential[3];

    double k_vec[3];

    double threeD_velocity[3] = {planet->velocity[0], planet->velocity[1], 0};

    crossProduct(position_vec, threeD_velocity, k_vec);
    if(k_vec[2] > 0){
        k_vec[2] = 1;
    }else if(k_vec[2] < 0){
        k_vec[2] = -1;
    }
    
    crossProduct(k_vec, position_vec, tangential); //two unit vectors crossed is still normalized 

    double twoD_tangential[2] = { tangential[0], tangential[1] };
    double tangential_velocity = dotproduct(twoD_tangential, planet->velocity, 2);

    double delta = k_vec[2] * tangential_velocity * TIME_STEP / planet->radius; //i.e. arc / radius = angle 
    planet->angle = planet->angle + delta;
    
    update_position_given_polar(planet->position, planet->angle, planet->radius);
    update_velocity(planet, planet->angle, sun_mass);

}


void init_planet(struct Kepler_Variables* planet, int planetColor, double mass, double* position, double* velocity, double* sun_pos) {

    double defaultpos[2] = { 1, 0 };

    planet->mass = mass;
    planet->position[0] = position[0];
    planet->position[1] = position[1];
    planet->velocity[0] = velocity[0];
    planet->velocity[1] = velocity[1];
    planet->angle = angle_between_2D_vectors(position, defaultpos);
    planet->previous_position[0] = 0;
    planet->previous_position[1] = 0;
    planet->active = true;
    planet->planetColor = planetColor;
    full_Kepler(planet, sun_pos, SUN_MASS);
}


//draw////////////////////////////////////////////////


void clearPlanet(double, double, uint16_t, int);
void drawPlanet(double x_d, double y_d, const uint16_t colour[BOXPIXELS][BOXPIXELS], int size, int rotation);
void redoBackground(double x_d, double y_d, const uint16_t bg[YRES][XRES], int size);


void drawPlanet(double x_d, double y_d, const uint16_t colour[BOXPIXELS][BOXPIXELS], int size, int rotation) {

    int x = (int)x_d;
    int y = (int)y_d;

    for (int i = -size / 2; i <= size/2; i++) {
        for (int j = -size / 2; j <= size/2; j++) {
            if ((i) * (i)+(j) * (j) <= size / 2 * size / 2)
                if (i + x < XRES && i + x > 0 && j + y > 0 && j + y < YRES)
                    plot_pixel(i + x, j + y, colour[(j + size / 2 ) % BOXPIXELS][i + size / 2 + rotation]); //horizontal rotation modulus so that index = BOXPIXELS becomes the 0th
        }
    }
}

void drawStar(double x_d, double y_d, const uint16_t colour[BOXPIXELS*2][BOXPIXELS*2], int size, int rotation) {

    int x = (int)x_d;
    int y = (int)y_d;

    for (int i = -size / 2; i <= size/2; i++) {
        for (int j = -size / 2; j <= size/2; j++) {
            if ((i) * (i)+(j) * (j) <= size / 2 * size / 2)
                if (i + x < XRES && i + x > 0 && j + y > 0 && j + y < YRES)
                    plot_pixel(i + x, j + y, colour[(j + size / 2 ) % BOXPIXELS * 2][i + size / 2 + rotation]); //horizontal rotation modulus so that index = BOXPIXELS becomes the 0th
        }
    }
}

void clearPlanet(double x_d, double y_d, uint16_t colour, int size) {

    int x = (int)x_d;
    int y = (int)y_d;

    for (int i = -size / 2; i <= size/2; i++) {
        for (int j = -size / 2; j <= size/2; j++) {
            if ((i) * (i)+(j) * (j) <= size / 2 * size / 2)
                if (i + x < XRES && i + x > 0 && j + y > 0 && j + y < YRES)
                    plot_pixel(i + x, j + y, colour); //horizontal rotation modulus so that index = BOXPIXELS becomes the 0th
        }
    }
}

void redoBackground(double x_d, double y_d, const uint16_t bg[YRES][XRES], int size) {
    int x = (int)x_d;
    int y = (int)y_d;

    for (int i = -size / 2; i <= size/2; i++) {
        for (int j = -size / 2; j <= size/2; j++) {
            if ((i) * (i)+(j) * (j) <= size / 2 * size / 2)
                if (i + x < XRES && i + x > 0 && j + y > 0 && j + y < YRES)
                    plot_pixel(i + x, j + y, bg[j + y ][i + x ]);
        }
    }
}


///////////////END OF DRAW///////////////////////

//integration functs

void user_selects_planet(Kepler_Variables* planet_ptr[20], int planetColor, double mass, double positionX, double positionY, double velocityX, double velocityY, double sun_pos[2]) {
    planet_ptr[NUM_PLANETS] = (struct Kepler_Variables*)malloc(sizeof(Kepler_Variables));
    double position[2] = { positionX, positionY };
    double velocity[2] = { velocityX, velocityY };
    init_planet(planet_ptr[NUM_PLANETS], planetColor, mass, position, velocity, sun_pos);
    NUM_PLANETS++;
}

void free_planets(Kepler_Variables** planetArr) {
    for (int i = 0; i < NUM_PLANETS; i++) {
        free(planetArr[i]);
    }
}


int main(void){

	short int colour = 0x39ff; //mouse colinstan\

	// // vga initialization
    *(pixel_ctrl_ptr + 1) = (int) &Buffer1; 
    wait_for_vsync();
    pixel_buffer_start = *pixel_ctrl_ptr;
    clear_screen(); 
    *(pixel_ctrl_ptr + 1) = (int) &Buffer2;
    pixel_buffer_start = *(pixel_ctrl_ptr + 1); 
    clear_screen(); 
	
	char welcometext[] = "WELCOME TO SPACE . . . \0";
	char chooseplanet[] = "SELECT YOUR PLANET \0";
	char planetmass[] = "Planet mass [kg]: \0";
	char velocity_x[] = "X Velocity [m/s]: \0";
	char velocity_y[] = "Y Velocity [m/s]: \0";
	char continuemsg[] = "PROCEED >> \0";
	char seldistance[] = "CLICK TO PLACE PLANET \0";

	char key0[8] = "[KEY0]\0";
	char key1[8] = "[KEY1]\0";
	char key2[8] = "[KEY2]\0";
	char key3[8] = "[KEY3]\0";
    start_screen = true;

    // clear all char buffers from last compile
    video_text_instant(29,29,clear);
    video_text_instant(31,5,clear);
    video_text_instant(21, 30, clear);
	video_text_instant(53, 30, clear);
	video_text_instant(21, 55, clear);
	video_text_instant(53, 55, clear);
	video_text_instant(31, 40, clear);
	video_text_instant(31, 50, clear);
	video_text_instant(31, 32, clear);
	video_text_instant(39, 35, clear);
	video_text_instant(39, 44, clear);
	video_text_instant(39, 53, clear);
	video_text_instant(69, 56, clear);
	video_text_instant(30, 1, clear);
	video_text_instant(30, 1, clear);

    double r_sun_pos[2] = { 0,0 };
    struct Kepler_Variables* planet_pointer_arr[20];

    // PS2 mouse initialization
    PS2_Mouse data = {0, 0, 0};

    PS2_send(0xFF);
    while(PS2_get() != 0xFA);
    
    while(PS2_get() != 0xAA);

    PS2_send(0xF0);
    while(PS2_get() != 0xFA);

    PS2_send(0xF4);
    while(PS2_get() != 0xFA);

    // drawPNG(startscreen);
    struct Point mouse = {320/2,240/2};
	mouse.prevx = 0;
	mouse.prevy = 0;

    while(1){

        short int colour = 0xFFFF;
        PS2_send(0xEB);
        while(PS2_get() != 0xFA);
        
        uint8_t status = PS2_get();
        int8_t x_del = PS2_get();
        int8_t y_del = PS2_get();

        christina:

		//delete mouse, replace with background
		if(!doneDraw2){//in start screen
			redraw_background_for_mouse(mouse, startingscreen);
		}
		else if(!doneDraw3){ // in select screen
			redraw_background_for_mouse(mouse, selectplanet);
		}
		else if(onblueplanet){
			redraw_background_for_mouse(mouse, blueplanetscreen);
		}
		else if(oniceplanet){
			redraw_background_for_mouse(mouse, iceplanetscreen);
		}
		else if(onredplanet){
			redraw_background_for_mouse(mouse, redplanetscreen);
		}
		else if(ongreenplanet){
			redraw_background_for_mouse(mouse, greenplanetscreen);
		}
		else if(onselectdistance){
			redraw_background_for_mouse(mouse, selectdistance);
		}

    	//old = new
		mouse.prevx = mouse.x;
		mouse.prevy = mouse.y;
	   
		//increment
        mouse.x += x_del;
        mouse.y -= y_del;


        HEX_PS2(status, mouse.x, mouse.y);

        if(status == 0x9){
            mouseClicked = true;
            clickPoint.x = mouse.x;
            clickPoint.y = mouse.y;
        }
        if(status != 0x9 && mouseClicked){
            negEdgeClicked = true;
        }


        if(start_screen){
            if(!doneDraw1){
                drawPNG(startingscreen);
				swapBuffer();
				drawPNG(startingscreen);
				swapBuffer();

                doneDraw1 = true;
            }
            if(negEdgeClicked){
                start_screen = false;
                welcome_screen = true;
                mouseClicked = false;
                negEdgeClicked = false;
            }
        }

        else if(welcome_screen){
            drawPNG(galaxy);
			swapBuffer();
			drawPNG(galaxy);
			swapBuffer();

            video_text_delay(29, 29, welcometext);
            video_text_delay(29, 29, clear);
            welcome_screen = false;
            select_screen = true;
        }

        else if(select_screen){
            if(!doneDraw2){
                drawPNG(selectplanet);
				swapBuffer();
				drawPNG(selectplanet);
				swapBuffer();
                video_text_delay(31, 5, chooseplanet);
                doneDraw2 = true;
            }
            if(negEdgeClicked){
                
                if((clickPoint.x > 63 && clickPoint.x < 128) && (clickPoint.y > 44 && clickPoint.y < 110)){
					onredplanet = false;
					ongreenplanet = false;
					oniceplanet = false;
					onblueplanet = false;
					video_text_instant(31,5,clear);
                    drawPNG(blueplanetscreen);
					swapBuffer();
					drawPNG(blueplanetscreen);
                    negEdgeClicked = false;
                    mouseClicked = false;
					doneDraw3 = true;
					onblueplanet = true;
					customization_screen = true;
                    select_screen = false;
					//redraw_background_for_mouse(mouse, blueplanetscreen);
                    
                }
				else if((clickPoint.x > 191 && clickPoint.x < 255) && (clickPoint.y > 44 && clickPoint.y < 110)){
					onredplanet = false;
					ongreenplanet = false;
					oniceplanet = false;
					onblueplanet = false;
					video_text_instant(31,5,clear);
                    drawPNG(iceplanetscreen);
					swapBuffer();
					drawPNG(iceplanetscreen);
                    negEdgeClicked = false;
                    mouseClicked = false;
					doneDraw3 = true;
					oniceplanet = true;
					customization_screen = true;
                    select_screen = false;
					//redraw_background_for_mouse(mouse, iceplanetscreen);
                    
				}
				else if((clickPoint.x > 63 && clickPoint.x < 128) && (clickPoint.y > 143 && clickPoint.y < 210)){
					onredplanet = false;
					ongreenplanet = false;
					oniceplanet = false;
					onblueplanet = false;
					video_text_instant(31,5,clear);
                    drawPNG(redplanetscreen);
					swapBuffer();
					drawPNG(redplanetscreen);
                    negEdgeClicked = false;
                    mouseClicked = false;
					doneDraw3 = true;
					onredplanet = true;
					customization_screen = true;
                    select_screen = false;
					//redraw_background_for_mouse(mouse, redplanetscreen);
                    
				}
				else if((clickPoint.x > 191 && clickPoint.x < 255) && (clickPoint.y > 143 && clickPoint.y < 210)){
					onredplanet = false;
					ongreenplanet = false;
					oniceplanet = false;
					onblueplanet = false;
					video_text_instant(31,5,clear);
                    drawPNG(greenplanetscreen);
					swapBuffer();
					drawPNG(greenplanetscreen);
                    negEdgeClicked = false;
                    mouseClicked = false;
					doneDraw3 = true;
					ongreenplanet = true;
					customization_screen = true;
                    select_screen = false;
					//redraw_background_for_mouse(mouse, greenplanetscreen);
                    
				}
            }
        }
		else if(customization_screen){
			int planetMass;
			int xVel;
			int yVel;

			if(!doneDraw4){
				planetIndex++;
				planets[planetIndex] = (planetPassIn){1000, 0, 10000, 0,0,0,0};

				planetMass = planets[planetIndex].mass;
				char stringMass[25];
				itoa(planetMass, stringMass, 10);
                strcat(stringMass, "\0");
                

				xVel = planets[planetIndex].x_velocity;
				char stringXVel[25];
				itoa(xVel, stringXVel, 10);
                strcat(stringXVel, "\0");


				yVel = planets[planetIndex].y_velocity;
				char stringYVel[25];
				itoa(yVel, stringYVel, 10);
                strcat(stringYVel, "\0");

				video_text_instant(31, 32, planetmass);
				video_text_instant(31, 40, velocity_x);
				video_text_instant(31, 50, velocity_y);

				video_text_instant(39, 35, stringMass);
				video_text_instant(39, 44, stringXVel);
				video_text_instant(39, 53, stringYVel);
				video_text_instant(69, 56, continuemsg);

				doneDraw4 = true;
			}
			// handle arrow clicking
			if(negEdgeClicked){
				// top left arrow mass
				if((clickPoint.x > 98 && clickPoint.x < 109) && (clickPoint.y > 139 && clickPoint.y < 144)){
					planets[planetIndex].mass -= 100;
					char stringMass[25];
					itoa((int)planets[planetIndex].mass, stringMass, 10);
                    strcat(stringMass, "\0");
                    video_text_instant(39, 35, "    ");
					video_text_instant(39, 35, stringMass);
                    
				}
				// top right
				else if((clickPoint.x > 201 && clickPoint.x < 214) && (clickPoint.y > 139 && clickPoint.y < 144)){
					planets[planetIndex].mass += 100;
					char stringMass[25];
					itoa((int)planets[planetIndex].mass, stringMass, 10);
                    strcat(stringMass, "\0");
                    video_text_instant(39, 35, "    ");
					video_text_instant(39, 35, stringMass);
				}
				// middle left
				else if((clickPoint.x > 98 && clickPoint.x < 109) && (clickPoint.y > 177 && clickPoint.y < 183)){
					planets[planetIndex].x_velocity -= 500;
					char stringxvel[25];
					itoa((int)planets[planetIndex].x_velocity, stringxvel, 10);
                    video_text_instant(39, 44, "     ");
					video_text_instant(39, 44, stringxvel);
				}
				// middle right
				else if((clickPoint.x > 201 && clickPoint.x < 214) && (clickPoint.y > 177 && clickPoint.y < 183)){
					planets[planetIndex].x_velocity += 500;
					char stringxvel[25];
					itoa((int)planets[planetIndex].x_velocity, stringxvel, 10);
                    video_text_instant(39, 44, "     ");
					video_text_instant(39, 44, stringxvel);
				}
				// bottom left
				else if((clickPoint.x > 98 && clickPoint.x < 109) && (clickPoint.y > 213 && clickPoint.y < 220)){
					planets[planetIndex].y_velocity -= 500;
					char *stringYVel;
					itoa( (int) planets[planetIndex].y_velocity, stringYVel, 10);
                    video_text_instant(39, 53, "     ");
					video_text_instant(39, 53, stringYVel);
				}
				// bottom right
				else if((clickPoint.x > 201 && clickPoint.x < 214) && (clickPoint.y > 213 && clickPoint.y < 220)){
					planets[planetIndex].y_velocity += 500;
					char *stringYVel;
					itoa( (int)planets[planetIndex].y_velocity, stringYVel, 10);
                    video_text_instant(39, 53, "     ");
					video_text_instant(39, 53, stringYVel);
				}
				// HANDLE PROCEED BUTTON
				else if((clickPoint.x > 276 && clickPoint.x < 314) && (clickPoint.y > 223 && clickPoint.y < 228)){
					video_text_instant(31, 32, clear);
					video_text_instant(31, 40, clear);
					video_text_instant(31, 50, clear);
					video_text_instant(39, 35, clear);
					video_text_instant(39, 44, clear);
					video_text_instant(39, 53, clear);
					video_text_instant(69, 56, clear);
					onselectdistance = true;
					customization_screen = false;
					select_distance_screen = true;
				}
				

				printf("x:");
            	printf("%d",mouse.x);
            	printf("\n");
           		printf("y:");
            	printf("%d", mouse.y);
            	printf("\n");
            	printf("clickstatus:");
            	printf("%d", status);
            	printf("\n");
            	printf("--");
            	printf("\n");
				negEdgeClicked = false;
            	mouseClicked = false;
			}
		}
		else if(select_distance_screen){

			if(!doneDraw5){
                drawPNG(selectdistance);
				swapBuffer();
				drawPNG(selectdistance);
				video_text_delay(30, 1, seldistance);
				doneDraw5 = true;
				if(oniceplanet){
					planets[planetIndex].planetColor = 0;
				}
				else if(onredplanet){
					planets[planetIndex].planetColor = 1;
				}
				else if(ongreenplanet){
					planets[planetIndex].planetColor = 2;
				}
				else if(onblueplanet){
					planets[planetIndex].planetColor = 3;
				}
                ongreenplanet = false;
                oniceplanet = false;
                onblueplanet = false;
                onredplanet = false;
			}

			if(negEdgeClicked){
                video_text_instant(30, 1, clear);
				int clickx = clickPoint.x - 160; 
				int clicky = clickPoint.y - 120;
                planets[planetIndex].x_pos = clickx;
                planets[planetIndex].y_pos = clicky;
				planets[planetIndex].distanceFromSun = sqrt((clickx*clickx) + (clicky*clicky));
				printf("%f", planets[planetIndex].distanceFromSun);
				negEdgeClicked = false;
            	mouseClicked = false;
                select_distance_screen = false;
                break;
			}

		}
		//draw mouse
		draw_point(mouse, colour);
        swapBuffer();
    }
        

    //init planet and sun////////////////////////////////
    

    printf("===========hi=======\n");
    printf("%f \n", planets[planetIndex].mass);
    printf("%f \n", planets[planetIndex].x_pos);
    printf("%f \n", planets[planetIndex].y_pos);
    printf("%f \n", planets[planetIndex].x_velocity);
    printf("%f \n", planets[planetIndex].y_velocity);
    user_selects_planet(planet_pointer_arr, planets[planetIndex].planetColor, planets[planetIndex].mass / 1000, planets[planetIndex].x_pos, planets[planetIndex].y_pos, planets[planetIndex].x_velocity/100000, planets[planetIndex].y_velocity/100000, r_sun_pos);
    

    int rotation = 0;
    bool showTrace = true;


    //PLANET MOVING SUPERLOOP
    while (simulation_time < 400000)
    {
        
        PS2_send(0xEB);
        while(PS2_get() != 0xFA);
        
        uint8_t status = PS2_get();
        
        HEX_PS2(status, mouse.x, mouse.y);

        if(status == 0x9){
            mouseClicked = true;
            clickPoint.x = mouse.x;
            clickPoint.y = mouse.y;
        }
        if(status != 0x9 && mouseClicked){
            negEdgeClicked = true;
        }

        
       
        drawStar(160, 120, sunsquare, BOXPIXELS * 2, rotation);

        //clear planets
        for (int i = 0; i < NUM_PLANETS; i++) {
            if (simulation_time != 0)
            {
                redoBackground(planet_pointer_arr[i]->previous_position[0] + XRES/2, planet_pointer_arr[i]->previous_position[1] + YRES / 2, galaxy, BOXPIXELS);

                if (showTrace) {
                    plot_pixel((int) (planet_pointer_arr[i]->previous_position[0] + XRES / 2 - planet_pointer_arr[i]->velocity[0] * 30),  (int) (planet_pointer_arr[i]->previous_position[1] + YRES / 2 - planet_pointer_arr[i]->velocity[1] * 30), 0xFFFF);
                }
            }
        }

        //save prev positions
        for (int i = 0; i < NUM_PLANETS; i++) {
            if (planet_pointer_arr[i]->active)
            {
                planet_pointer_arr[i]->previous_position[0] = planet_pointer_arr[i]->position[0];
                planet_pointer_arr[i]->previous_position[1] = planet_pointer_arr[i]->position[1];

            }
        }

        //increment rotation
        if (rotation == BOXPIXELS) {
            rotation = 0;
        }
        else {
            rotation++;
        }

        for (int i = 0; i < NUM_PLANETS; i++) {
            //inc position
            if (planet_pointer_arr[i]->active) {
                orbit(planet_pointer_arr[i], SUN_MASS);

                if (planet_pointer_arr[i]->position[0] + XRES / 2 - (int)(BOXPIXELS / 2) >= XRES || planet_pointer_arr[i]->position[1] + YRES / 2 - (int)(BOXPIXELS / 2) >= YRES) {
                    planet_pointer_arr[i]->active = false;
                    printf("planet is out of bounds =======================\n");
                    clearPlanet(planet_pointer_arr[i]->position[0] + XRES / 2, planet_pointer_arr[i]->position[1] + YRES / 2, 0x0, BOXPIXELS);

                }

                if (planet_pointer_arr[i]->radius < 0) {
                    printf("negative radius error");
                    break;
                }
            }
        }


        //redraw
        for (int i = 0; i < NUM_PLANETS; i++) {
            //draw
            if (planet_pointer_arr[i]->active) {
                    switch(planet_pointer_arr[i]->planetColor){
                        case 0:
                        drawPlanet(planet_pointer_arr[i]->position[0] + XRES / 2, planet_pointer_arr[i]->position[1] + YRES / 2, smalliceplanet, BOXPIXELS, rotation);
                        break;
                        case 1:
                        drawPlanet(planet_pointer_arr[i]->position[0] + XRES / 2, planet_pointer_arr[i]->position[1] + YRES / 2, smallredplanet, BOXPIXELS, rotation);
                        break;
                        case 2:
                        drawPlanet(planet_pointer_arr[i]->position[0] + XRES / 2, planet_pointer_arr[i]->position[1] + YRES / 2, smallgreenplanet, BOXPIXELS, rotation);
                        break;
                        case 3:
                        drawPlanet(planet_pointer_arr[i]->position[0] + XRES / 2, planet_pointer_arr[i]->position[1] + YRES / 2, smallblueplanet, BOXPIXELS, rotation);
                        break;
                        default:
                        drawPlanet(planet_pointer_arr[i]->position[0] + XRES / 2, planet_pointer_arr[i]->position[1] + YRES / 2, smalliceplanet, BOXPIXELS, rotation);
                    }
                //printf("position x: %d  position y: %d  time: %f  radius: %d  angle: %d  ========== \n", (int)(planet_pointer_arr[i]->position[0] / DOWNSCALE), (int)(planet_pointer_arr[i]->position[1] / DOWNSCALE), simulation_time, (int)planet_pointer_arr[i]->radius, (int)rad_to_angle(planet_pointer_arr[i]->angle));
            }
        }

        //check for cancel
        if(negEdgeClicked){
                mouse.x = XRES/2;
                mouse.y = YRES/2;
                mouseClicked = false;
                negEdgeClicked = false;
                doneDraw2 = false;
                doneDraw3 = false;
                doneDraw4 = false;
                doneDraw5 = false;
                select_screen = true;

				goto christina;
			}

        swapBuffer();
        simulation_time += 10;

    }
    for (int i = 0; i < NUM_PLANETS; i++) {
        free(planet_pointer_arr[i]);
    }

        
        return 0;
    }


void swapBuffer(){
	wait_for_vsync();
    pixel_buffer_start = *(pixel_ctrl_ptr + 1); 
}


void PS2_send(uint8_t c){
    ps2->data = c;
    // return PS2_get() == 0xFA ? 0 : 1;
}

uint8_t PS2_get(){
    uint32_t data, rvalid;
    do{
        data = ps2->data;
        rvalid = (data >> 15) & 1;
    }while(!rvalid);
    return data & 0xFF;
}

void HEX_PS2(char b1, char b2, char b3) {
    volatile int * HEX3_HEX0_ptr = (int *)HEX3_HEX0_BASE;
    volatile int * HEX5_HEX4_ptr = (int *)HEX5_HEX4_BASE;
    /* SEVEN_SEGMENT_DECODE_TABLE gives the on/off settings for all segments in
    * a single 7-seg display in the DE1-SoC Computer, for the hex digits 0 - F
    */
    unsigned char seven_seg_decode_table[] = {
        0x3F, 0x06, 0x5B, 0x4F, 0x66, 0x6D, 0x7D, 0x07,
        0x7F, 0x67, 0x77, 0x7C, 0x39, 0x5E, 0x79, 0x71};
    unsigned char hex_segs[] = {0, 0, 0, 0, 0, 0, 0, 0};
    unsigned int shift_buffer, nibble;
    unsigned char code;
    int i;
    shift_buffer = (b1 << 16) | (b2 << 8) | b3;
    for (i = 0; i < 6; ++i) {
        nibble = shift_buffer & 0x0000000F; // character is in rightmost nibble
        code = seven_seg_decode_table[nibble];
        hex_segs[i] = code;
        shift_buffer = shift_buffer >> 4;
    }
    /* drive the hex displays */
    *(HEX3_HEX0_ptr) = *(int *)(hex_segs);
    *(HEX5_HEX4_ptr) = *(int *)(hex_segs + 4);
}


void wait_for_vsync(){
    *pixel_ctrl_ptr = 1;
    while ((*pixel_status & 1) != 0);
}


struct timer_t * const timer = (struct timer_t *) 0xFF202000;

void waitasec(int pow_fraction) {
	   unsigned int t = TIMERSEC >> pow_fraction;
       timer->control = 0x8; // stop the timer
       timer->status = 0; // reset TO
       timer->periodlo = (t & 0x0000FFFF);
       timer->periodhi = (t & 0xFFFF0000) >> 16;
       timer->control = 0x4;
       while ((timer->status & 0x1) == 0);
          timer->status = 0; // reset TO
}

void video_text_delay(int x, int y, char * text_ptr) {
	int offset;

	volatile char * character_buffer = (char *)FPGA_CHAR_BASE; // video character buffer
	/* assume that the text string fits on one line */
	offset = (y << 7) + x;
	while (*(text_ptr)) {
		*(character_buffer + offset) = *(text_ptr); // write to the character buffer
		++text_ptr;
		++offset;
		waitasec(1);
	}
}

void video_text_instant(int x, int y, char * text_ptr) {
	int offset;
	volatile char * character_buffer = (char *)FPGA_CHAR_BASE; // video character buffer
	/* assume that the text string fits on one line */

	offset = (y << 7) + x;
	while (*(text_ptr)) {
		*(character_buffer + offset) = *(text_ptr); // write to the character buffer
		++text_ptr;
		++offset;
	}
}



void drawPNG(const uint16_t arr[YRES][XRES]){
    for(int x = 0; x < XRES; x++){
        for(int y = 0; y < YRES; y++){
            plot_pixel(x,y,arr[y][x]);
        }
    }
}


void redraw_background_for_mouse(struct Point p , const uint16_t arr[YRES][XRES]){
      int x_min = p.prevx - POINT_PADDING;
    int y_min = p.prevy - POINT_PADDING;
    int x_max = p.prevx + POINT_PADDING;
    int y_max = p.prevy + POINT_PADDING;

    x_min = 0 > x_min ? 0 : x_min;
    y_min = 0 > y_min ? 0 : y_min;
    x_max = 319 < x_max ? 319 : x_max;
    y_max = 239 < y_max ? 239 : y_max;

    for(int x = x_min; x <= x_max; x++)
        for(int y = y_min; y <= y_max; y++)
            if(x == p.prevx || y == p.prevy)
            plot_pixel(x,y,arr[y][x]);
}

void clear_screen(){
    int x, y;
	for (x = 0; x < 320; x++)
       for (y = 0; y < 240; y++) 
          plot_pixel(x, y, 0x00);
}

void drawBox( int x, int y, short colour){
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            plot_pixel(i + x, j + y, colour);
        }
    }
}

void plot_pixel(int x, int y, uint16_t line_color)
{

    volatile short int* one_pixel_address;

    if (x < XRES && x > 0 && y < YRES && y > 0) {
        one_pixel_address = pixel_buffer_start + (y << 10) + (x << 1);
        *one_pixel_address = line_color;
    }
}

void draw_point(struct Point p, short int colour){

    int x_min = p.x - POINT_PADDING;
    int y_min = p.y - POINT_PADDING;
    int x_max = p.x + POINT_PADDING;
    int y_max = p.y + POINT_PADDING;

    x_min = 0 > x_min ? 0 : x_min;
    y_min = 0 > y_min ? 0 : y_min;
    x_max = 319 < x_max ? 319 : x_max;
    y_max = 239 < y_max ? 239 : y_max;

    for(int x = x_min; x <= x_max; x++)
        for(int y = y_min; y <= y_max; y++)
        if(x == p.x || y == p.y)
            plot_pixel(x,y,colour);
}


void swap(int* x, int* y){
    int temp;
    temp = *x;
    *x = *y;
    *y = temp;
}
