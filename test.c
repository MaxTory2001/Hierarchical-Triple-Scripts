#include <stdio.h>

int main(int argc, char* argv[]) {
    printf("Received %d arguments\n", argc);
    int i = 0;
    for (i = 0; i < argc; i++) {
	printf("arg%1d = %s\n", i, argv[i]);
    }
}
