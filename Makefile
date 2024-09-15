# make file for swan
compiler = gfortran
target = main
src_dir = src
mod_dir = mod
obj_dir = obj
mods = helper.f90 timestep.f90
f90 = main.f90
srcs := $(mods) $(f90)
srcs := $(addprefix $(src_dir)/, $(srcs))
objs := $(patsubst $(src_dir)/%.f90,$(obj_dir)/%.o, $(srcs))


#define check_and_create_folder
#	@echo "hah"
#	FOLDER_NAME = $(1)
#	ifeq ($(wildcard $(FOLDER_NAME)),)
#		@echo "Creating $(FOLDER_NAME) folder"
#		@mkdir -p $(FOLDER_NAME)
#	else
#		@echo "$(FOLDER_NAME) folder already exists"
#	endif
#endef


$(target):create_folders $(objs)
	$(compiler) $(objs) -o $(target)

$(obj_dir)/%.o:$(src_dir)/%.f90
	$(compiler) -c $< -o $@ -J $(mod_dir)

create_folders:
ifeq ($(wildcard $(mod_dir)),)
	@echo "Creating $(mod_dir) folder"
	@mkdir -p $(mod_dir)
endif
ifeq ($(wildcard $(obj_dir)),)
	@echo "Creating $(obj_dir) folder"
	@mkdir -p $(obj_dir)
endif
#	$(call check_and_create_folder, $(mod_dir))
#	$(call check_and_create_folder, $(obj_dir))

clean:
	rm -rf $(mod_dir) $(obj_dir)

cleanoutput:
	rm -rf output

cleanall:
	rm -rf $(target) output $(mod_dir) $(obj_dir)