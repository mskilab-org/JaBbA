#!/bin/bash

install_cplex() {
    local CPLEX_BIN_NAME=$1
    local CPLEX_INSTALL_PATH=$2
    local CPLEX_DIR="$CPLEX_INSTALL_PATH/CPLEX_Studio"
    local CONTAINER_NAME="jabba_cplex_container"
    local NEW_IMAGE_NAME="jabba_cplex_image"
    local IMAGE_NAME="mskilab/jabba"

    echo "Installing CPLEX..."
    docker pull $IMAGE_NAME:latest
    CONTAINER_ID=$(docker run -itd --rm --platform linux/amd64 --name $CONTAINER_NAME $IMAGE_NAME:latest)
    docker exec $CONTAINER_ID mkdir -p /opt/cplex_studio/
    docker cp ./$CPLEX_BIN_NAME $CONTAINER_ID:/opt/cplex_studio/
    docker exec $CONTAINER_ID chmod 777 /opt/cplex_studio/$CPLEX_BIN_NAME
    echo "Installing CPLEX, this may take a few minutes..."
    docker exec $CONTAINER_ID /opt/cplex_studio/$CPLEX_BIN_NAME -i silent -DLICENSE_ACCEPTED=TRUE -DUSER_INSTALL_DIR=$CPLEX_INSTALL_PATH
    echo "CPLEX has been installed, setting CPLEX_DIR environment variable..."
    ENV_CHANGE="ENV CPLEX_DIR=$CPLEX_DIR"
    docker commit --change="$ENV_CHANGE" $CONTAINER_ID $NEW_IMAGE_NAME
    docker stop $CONTAINER_ID
    docker save $NEW_IMAGE_NAME > jabba_cplex.img
    echo "CPLEX has been installed, CPLEX_DIR set, and the new image has been saved as jabba_cplex.img"
}

install_gurobi() {
    local GUROBI_TAR_NAME=$1
    local GUROBI_INSTALL_PATH=$2
    local GUROBI_DIR="$GUROBI_INSTALL_PATH/linux64"
    local CONTAINER_NAME="jabba_gurobi_container"
    local NEW_IMAGE_NAME="jabba_gurobi_image"
    local IMAGE_NAME="mskilab/jabba"

    echo "Installing GUROBI..."
    docker pull $IMAGE_NAME:latest
    CONTAINER_ID=$(docker run -itd --rm --platform linux/amd64 --name $CONTAINER_NAME $IMAGE_NAME:latest)
    docker exec $CONTAINER_ID mkdir -p $GUROBI_INSTALL_PATH
    docker cp ./$GUROBI_TAR_NAME $CONTAINER_ID:$GUROBI_INSTALL_PATH
    echo "Installing GUROBI, this may take a few minutes..."
    docker exec $CONTAINER_ID tar -xvf $GUROBI_INSTALL_PATH/$GUROBI_TAR_NAME -C $GUROBI_INSTALL_PATH
    echo "GUROBI has been installed, setting GUROBI_HOME environment variable..."
    ENV_CHANGE="ENV GUROBI_HOME=$GUROBI_DIR"
    docker commit --change="$ENV_CHANGE" $CONTAINER_ID $NEW_IMAGE_NAME
    docker stop $CONTAINER_ID
    docker save $NEW_IMAGE_NAME > jabba_gurobi.img
    echo "GUROBI has been installed, GUROBI_HOME set, and the new image has been saved as jabba_gurobi.img"
}

# Main function to handle user input and call the appropriate installation function
main() {
    echo "Select the software to install:"
    echo "1. CPLEX"
    echo "2. GUROBI"
    read -p "Enter your choice (1 or 2): " choice

    case $choice in
        1)
            read -p "Enter the CPLEX binary name (e.g., install_cplex.bin): " cplex_bin_name
            read -p "Enter the path where CPLEX should be installed (e.g., /opt/cplex): " cplex_install_path
            install_cplex "$cplex_bin_name" "$cplex_install_path"
            ;;
        2)
            read -p "Enter the GUROBI tar.gz name (e.g., gurobi.tar.gz): " gurobi_tar_name
            read -p "Enter the path where GUROBI should be installed (e.g., /opt/gurobi): " gurobi_install_path
            install_gurobi "$gurobi_tar_name" "$gurobi_install_path"
            ;;
        *)
            echo "Invalid choice. Please enter 1 for CPLEX or 2 for GUROBI."
            ;;
    esac
}

# Call the main function
main
