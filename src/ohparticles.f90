#define OH_LIB_LEVEL 3

#ifndef OH_PBUF_SIZE
#define OH_PBUF_SIZE 16384
#endif
module m_ohparticles
    use ohhelp2, only: oh2_max_local_particles
    use oh_type, only: oh_particle
    implicit none

    private
    public t_OhParticles, new_OhParticles

    integer, parameter :: PARTICLE_BUFSIZE = OH_PBUF_SIZE

    type :: t_OhParticles
        type(oh_particle), allocatable :: pbuf(:)
        integer, allocatable :: pbase(:)
        integer(kind=8), private :: max_nparticles
        integer :: max_local_particles

        integer, allocatable :: total_local_particles(:, :) ! (nspec, 2)

        integer :: nspecies

    contains

        procedure :: allocate_pbuf => ohparticles_allocate_pbuf
        procedure :: start_index => ohparticles_start_index
        procedure :: end_index => ohparticles_end_index

    end type

contains

    function new_OhParticles(nspecies, &
                             max_nparticles, &
                             nprocs) result(self)
        integer, intent(in) :: nspecies
        integer(kind=8), intent(in) :: max_nparticles
        integer, intent(in) :: nprocs
        type(t_OhParticles) :: self

        self%nspecies = nspecies
        self%max_nparticles = max_nparticles

        allocate (self%pbase(3))
        allocate (self%total_local_particles(nspecies, 2))
    end function

    subroutine ohparticles_allocate_pbuf(self, loadbalance_tolerance_percentage)
        class(t_OhParticles), intent(inout) :: self
        integer, intent(in) :: loadbalance_tolerance_percentage

        self%max_local_particles = oh2_max_local_particles(self%max_nparticles, &
                                                           loadbalance_tolerance_percentage, &
                                                           PARTICLE_BUFSIZE)
        allocate (self%pbuf(self%max_local_particles))
    end subroutine

    function ohparticles_start_index(self, ispec, ps) result(ret)
        class(t_OhParticles), intent(in) :: self
        integer, intent(in) :: ispec
        integer, intent(in) :: ps
        integer :: ret

        ret = self%pbase(ps) + 1 + sum(self%total_local_particles(1:ispec - 1, ps))
    end function

    function ohparticles_end_index(self, ispec, ps) result(ret)
        class(t_OhParticles), intent(in) :: self
        integer, intent(in) :: ispec
        !> 1: primary subdomain, 2: secondary subdomain
        integer, intent(in) :: ps
        integer :: ret

        ret = self%pbase(ps) + sum(self%total_local_particles(1:ispec, ps))
    end function

end module
