module m_particle_mover_factory
    use m_particle_mover, only: t_ParticleMover
    use m_particle_mover_boris, only: t_ParticleMoverBoris, new_ParticleMoverBoris

    private
    public t_ParticleMoverFactory, new_ParticleMoverFactory

    type t_ParticleMoverFactory
    contains

        procedure :: create_particle_mover => factory_create_particle_mover
    end type

contains

    function new_ParticleMoverFactory() result(obj)
        type(t_ParticleMoverFactory) :: obj
    end function

    subroutine factory_create_particle_mover(self, name, particle_mover)
        class(t_ParticleMoverFactory), intent(in) :: self
        character(len=*), intent(in) :: name
        class(t_ParticleMover), allocatable, intent(out) :: particle_mover

        select case (name)
        case ('boris')
            allocate (particle_mover, source=new_ParticleMoverBoris())
        case default
            ! TODO: エラー出力?
        end select
    end subroutine

end module
