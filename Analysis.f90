program Analysis
    implicit none
    integer::type !どの解析を行うか

    write(*,*) "What Analysis?"
    write(*,*) "0: Plane_Calculation"
    write(*,*) "1: Wing Analysis"
    write(*,*) "2: Tail Analysis"
    write(*,*) "3: Fin Analysis"
    write(*,*) "4: Plane Analysis"
    write(*,*) "5: Polar Analysis"
    write(*,*) "6: NLFS6"
    write(*,*) "7: GA_NLFS6"
    read(*,*) type

    if(type==0) call Test_calculation()
    if(type==1) call Wing_analysis()
    if(type==2) call Polar_analysis()
    if(type==3) call Polar_analysis()
    if(type==4) call Polar_analysis()
    if(type==5) call Polar_analysis()
    if(type==6) call Test_NLFS()
    if(type==7) call GA_NLFS()

end program 

subroutine Test_calculation()
    implicit none !暗黙の変数宣言を無効にする
    !Plane_calc用
    real(8),dimension(0:2)::Altitude_angle !姿勢角 [Φ,θ,Ψ] [deg]
    real(8),dimension(0:2)::Ground_speed !機体軸における対地速度ベクトル [u,v,w] [m/s]
    real(8),dimension(0:2)::Angular_velocity !機体軸における角速度ベクトル [p,q,r] [deg/s]
    real(8),dimension(0:2)::gust !機体軸における突風ベクトル [ug,vg,wg] [m/s]
    real(8),dimension(0:2)::Input !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg] [δe,δh,δr]
    real(8)::hE !高度[m]
    real(8)::Lift !揚力 [N]
    real(8)::Drag !抗力 [N]
    real(8),dimension(0:2)::Force !機体にはたらく外力 [X,Y,Z] [N]
    real(8),dimension(0:2)::Moment !機体にはたらくモーメント [L,M,N] [N*m]
    real(8)::epsilon
    real(8),dimension(0:2)::Acceleration !機体軸における対地加速度ベクトル [du/dt,dv/dt,dw/dt] [m/s]
    real(8),dimension(0:2)::Angular_acceleration !機体軸における角加速度ベクトル [dp/dt,dq/dt,dr/dt] [deg/s]
    !その他
    real(8)::Velocity
    real(8)::X,Y,Z,L,M,N
    real(8)::phi,theta,psi
    real(8)::de,dh,dr
    real(8)::p,q,r
    real(8)::u,v,w,ug,vg,wg
    real(8)::alpha,beta,gamma
    real(8)::pi,deg,rad !円周率，degへの変換，radへの変換
    real(8)::tmp
    !計算時間計測
    real(8)::CPU_time_s
    integer::initial_time,finish_time,time_rate,CPU_time_min !initial time, finish time, time rate

    !値の入力
    write(*,*) "Velocity?"
    read(*,*) Velocity
    write(*,*) "alpha?"
    read(*,*) alpha
    write(*,*) "beta?"
    read(*,*) beta

    call system_clock(initial_time) !開始時間の読み込み
    
    pi=4.0D0*atan(1.0D0) !円周率
    deg=(180.0D0/pi) !degへの変換
    rad=(pi/180.0D0) !radへの変換

    !配列の初期化
    Lift=0.0D0 !揚力 [N]
    Drag=0.0D0 !抗力 [N]
    Force=0.0D0 !機体にはたらく外力 [X,Y,Z] [N]
    Moment=0.0D0 !機体にはたらくモーメント [L,M,N] [N*m]
    Altitude_angle=0.0D0
    Ground_speed=0.0D0 !機体軸における対地速度ベクトル [m/s]
    Angular_velocity=0.0D0 !角速度ベクトル [deg/s]
    Input=0.0D0 !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg]
    gust=0.0D0 !突風ベクトル [m/s]

    !初期値の入力
    gust(0)=0.0D0 !ug [m/s]
    gust(1)=0.0D0 !vg [m/s]
    gust(2)=0.0D0 !wg [m/s]
    Ground_speed(0)=Velocity*cos(rad*alpha)*cos(rad*beta) !u [m/s]
    Ground_speed(1)=Velocity*sin(rad*beta) !v [m/s]
    Ground_speed(2)=velocity*sin(rad*alpha)*cos(rad*beta) !w [m/s]
    Angular_velocity(0)=0.000D0 !p [deg/s]
    Angular_velocity(1)=0.000D0 !q [deg/s]
    Angular_velocity(2)=0.000D0 !r [deg/s]
    Altitude_angle(0)=0.000D0 !phi [deg]
    Altitude_angle(1)=0.000D0 !theta [deg]
    Altitude_angle(2)=0.000D0 !psi [deg]
    Input(0)=0.0D0 !de [deg]
    Input(1)=0.000D0 !dh [-]
    Input(2)=0.0D0 !dr [deg]
    hE=10.0D0
    epsilon=0.0D0
    !call write1D(3,Altitude_angle)

    call Plane_calculation &
        (Altitude_angle,Ground_speed,Angular_velocity,gust,Input &
        ,hE,Lift,Drag,Force,Moment,epsilon,Acceleration,Angular_acceleration)
    
    !わかりやすいように変数変換
    u=Ground_speed(0);          w=Ground_speed(2);              v=Ground_speed(1)           !u,v,w [m/s]
    ug=gust(0);                 wg=gust(2);                     vg=gust(1)                  !u,v,w [m/s]
    phi=rad*Altitude_angle(0);  theta=rad*altitude_angle(1);    psi=rad*altitude_angle(2)   !Φ，θ，ψ [rad]を格納
    p=Angular_velocity(0);      q=Angular_velocity(1);          r=Angular_velocity(2)       !p,q,r [deg]       
    de=Input(0);                dh=Input(1);                    dr=Input(2)                 !操舵を格納
    
    if(1==1) then
        write(*,*)
        write(*,'(A12,F12.3,A12,F12.3)') "Velosity",velocity,"he",he
        write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "alpha",alpha,"beta",beta,"gamma",gamma
        write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "Lift",Lift,"Drag",Drag,"L/D",Lift/Drag
        write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "u",u,"v",v,"w",w
        write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "phi",deg*phi,"theta",deg*theta,"psi",deg*psi
        write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "p",p,"q",q,"r",r
        write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "de",de,"dh",dh,"dr",dr
        write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "X",Force(0),"Y",Force(1),"Z",Force(2)
        write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "L",Moment(0),"M",Moment(1),"N",Moment(2)
        write(*,'(A12,F12.5,A12,F12.5,A12,F12.5)')  &
        "du/dt",Acceleration(0),"dv/dt",Acceleration(1),"dw/dt",Acceleration(2)
        write(*,'(A12,F12.5,A12,F12.5,A12,F12.5)')  &
        "dp/dt",Angular_Acceleration(0),"dq/dt",Angular_Acceleration(1),"dr/dt",Angular_Acceleration(2)
        write(*,*)
    end if

    call system_clock(finish_time,time_rate) !終了時間の読み込み
    CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
    if(CPU_time_s<60) then
        write(*,'(A24,f12.3,A12)') "Computation Time",CPU_time_s,"[s]"
    else 
        CPU_time_min=int(CPU_time_s/60.0D0)
        CPU_time_s=CPU_time_s-60.D0*CPU_time_min
        write(*,'(A24,I12,A12,F12.3,A12)') "Computation Time",CPU_time_min,"[min]",CPU_time_s,"[s]"
    end if

end subroutine

subroutine Wing_analysis()
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan,sqrt,sin,cos,tan !プログラム中で使用する関数の名前を宣言
    !パラメータ
    real(8),parameter::g =9.795D0 !重力加速度 [m/s^2]
    real(8),parameter::alpha_min =0.0D0 !最小迎角 [deg]
    real(8),parameter::alpha_max =30.0D0 !最大迎角 [deg]
    real(8),parameter::da =0.5D0 !迎角刻み [deg]
    real(8),parameter::S=18.332D0 !主翼面積．あくまで概算なので適当な値を [m^2]
    real(8),parameter::error=0.00049D0 !許容誤差
    !Plane_calculationの引数
    real(8),dimension(0:2)::Altitude_angle !姿勢角 [Φ,θ,Ψ] [deg]
    real(8),dimension(0:2)::Ground_speed !機体軸における対地速度ベクトル [u,v,w] [m/s]
    real(8),dimension(0:2)::Angular_velocity !機体軸における角速度ベクトル [p,q,r] [deg/s]
    real(8),dimension(0:2)::Input !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg] [δe,δh,δr]
    real(8),dimension(0:2)::gust !機体軸における突風ベクトル [ug,vg,wg] [m/s]
    real(8),dimension(0:2)::Acceleration !機体軸における対地加速度ベクトル [du/dt,dv/dt,dw/dt] [m/s]
    real(8),dimension(0:2)::Angular_acceleration !機体軸における角加速度ベクトル [dp/dt,dq/dt,dr/dt] [deg/s]
    real(8)::hE !高度 [m]
    real(8)::Lift,Drag !揚力，抗力 [N]
    real(8),dimension(0:2)::Force !機体にはたらく外力 [X,Y,Z] [N]
    real(8),dimension(0:2)::Moment !機体にはたらくモーメント [L,M,N] [N*m]
    real(8)::epsilon !主翼位置での吹きおろし角 [deg]
    !機体データ
    integer::span_div,chord_div !スパン方向分割数
    real(8)::Mass !機体重量 [kg]
    real(8)::Ixx,Iyy,Izz,Ixz !慣性モーメント [kg*m^4]
    real(8)::CD_fuselage,S_fuselage,chord_fuselage !カウルの抗力係数,面積，コード長
    !その他
    integer::pilot !操縦方法(0:重心移動，1:エレベーター)
    integer::parameter !解析で変化させるパラメーター(0:alpha, 1:beta)
    real(8)::Drag_fuselage
    real(8)::CL,CD,Cm,CL_old
    real(8)::rho,mu,Re
    real(8)::Velocity
    real(8)::X,Y,Z,L,M,N
    real(8)::phi,theta,psi
    real(8)::de,dh,dr
    real(8)::p,q,r
    real(8)::u,v,w,ug,vg,wg
    real(8)::alpha,beta,gamma
    real(8)::chord_mac,at !MAC,水平尾翼揚力傾斜(オールフライング，対称翼を仮定) [N/deg]
    real(8)::lt,it,Lift_tail,M0_wing,Lift_wing,hac,hcg
    real(8)::pi,deg,rad !円周率，degへの変換，radへの変換
    real(8)::tmp
    !カウンター
    integer::i,j,k

    pi=4.0D0*atan(1.0D0) !円周率
    deg=(180.0D0/pi) !degへの変換
    rad=(pi/180.0D0) !radへの変換

    Altitude_angle=0.0D0 !姿勢角 [Φ,θ,Ψ] [deg]
    Ground_speed=0.0D0 !機体軸における対地速度ベクトル [u,v,w] [m/s]
    Angular_velocity=0.0D0 !機体軸における角速度ベクトル [p,q,r] [deg/s]
    Input=0.0D0 !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg] [δe,δh,δr]
    gust=0.0D0 !機体軸における突風ベクトル [ug,vg,wg] [m/s]
    Acceleration=0.0D0 !機体軸における対地加速度ベクトル [du/dt,dv/dt,dw/dt] [m/s]
    Angular_acceleration=0.0D0 !機体軸における角加速度ベクトル [dp/dt,dq/dt,dr/dt] [deg/s]

    !わかりやすいように変数変換
    u=Ground_speed(0);          w=Ground_speed(2);              v=Ground_speed(1)           !u,v,w [m/s]
    ug=gust(0);                 wg=gust(2);                     vg=gust(1)                  !u,v,w [m/s]
    phi=Altitude_angle(0);      theta=Altitude_angle(1);            psi=Altitude_angle(2)   !Φ，θ，ψ [rad]を格納
    p=Angular_velocity(0);      q=Angular_velocity(1);          r=Angular_velocity(2)       !p,q,r [deg]       
    de=Input(0);                dh=Input(1);                    dr=Input(2)                 !操舵を格納

    write(*,*) "Velocity?"
    read(*,*) Velocity
    write(*,*) "hE?"
    read(*,*) hE
    write(*,*) "parameter? 0:alpha 1:beta"
    read(*,*) parameter

    !迎角α，横滑り角β [deg]に初期値を入れる
    if(parameter==0) then
        write(*,*) "beta?"
        read(*,*) beta
    else if(parameter==1) then
        write(*,*) "alpha?"
        read(*,*) alpha
    else
        stop
    end if

    write(*,*)
    write(*,'(A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10)') &
        "alpha","beta","CL","CDi","CDv","CD","Cy","Cl","Cm","Cn","Cni","Qinf","Xcp","L/D"


    do i=0,int((alpha_max-alpha_min)/da)
        
        if(parameter==0) alpha=alpha_min+da*real(i)
        if(parameter==1) beta=alpha_min+da*real(i)

        Lift=0.0D0 !揚力 [N]
        Drag=0.0D0 !抗力 [N]
        Force=0.0D0 !機体にはたらく外力 [X,Y,Z] [N]
        Moment=0.0D0 !機体にはたらくモーメント [L,M,N] [N*m]

        !主翼計算    
        open(10,file="Wing_data.txt") !ファイルを開く
        read(10,*) span_div !スパン方向分割数，コード方向分割数
        read(10,*) rho,mu !空気密度 [kg/m^3]，動粘性係数 [m^2/s]
        read(10,*) hac,tmp !空力中心位置
        read(10,*) tmp !翼素間隔 [m]
        read(10,*) chord_mac,hcg !平均空力翼弦長 [m]，重心位置 [-]
        close(10) !ファイルを閉じる
        call Wing_calculation(Lift,Drag,Force,Moment,Velocity,alpha,beta,p,r,hE,dh,epsilon,span_div)
        if(isnan(Lift)) write(*,'(A24)') "Wing"
        Lift_wing=Lift
        M0_wing=Moment(1)-Lift_wing*chord_mac*(hcg-dh-hac)

        !CL,CD,Cm,Gammaを計算
        CL=Lift/(0.5D0*rho*Velocity**2*S)
        CD=Drag/(0.5D0*rho*Velocity**2*S)
        Cm=Moment(1)/(0.5D0*rho*Velocity**2*S*chord_mac)

        write(*,'(F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6)') &
            alpha,beta,CL,0.0D0,0.0D0,CD,0.0D0,0.0D0,Cm,0.0D0,0.0D0,Velocity,0.0D0,Lift/Drag        

    end do

end subroutine

subroutine Polar_analysis()
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan,sqrt,sin,cos,tan !プログラム中で使用する関数の名前を宣言
    !パラメータ
    real(8),parameter::g =9.795D0 !重力加速度 [m/s^2]
    real(8),parameter::alpha_min =-3.0D0 !最小迎角 [deg]
    real(8),parameter::alpha_max =20.0D0 !最大迎角 [deg]
    real(8),parameter::da =0.1D0 !迎角刻み [deg]
    real(8),parameter::Velocity_max=14.0D0 !最大速度 [m/s]
    real(8),parameter::Velocity_min=7.5D0 !失速速度 [m/s]
    real(8),parameter::S=18.800D0 !主翼面積．あくまで概算なので適当な値を [m^2]
    real(8),parameter::error=0.00049D0 !許容誤差
    !Plane_calculationの引数
    real(8),dimension(0:2)::Altitude_angle !姿勢角 [Φ,θ,Ψ] [deg]
    real(8),dimension(0:2)::Ground_speed !機体軸における対地速度ベクトル [u,v,w] [m/s]
    real(8),dimension(0:2)::Angular_velocity !機体軸における角速度ベクトル [p,q,r] [deg/s]
    real(8),dimension(0:2)::Input !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg] [δe,δh,δr]
    real(8),dimension(0:2)::gust !機体軸における突風ベクトル [ug,vg,wg] [m/s]
    real(8),dimension(0:2)::Acceleration !機体軸における対地加速度ベクトル [du/dt,dv/dt,dw/dt] [m/s]
    real(8),dimension(0:2)::Angular_acceleration !機体軸における角加速度ベクトル [dp/dt,dq/dt,dr/dt] [deg/s]
    real(8)::hE !高度 [m]
    real(8)::Lift,Drag !揚力，抗力 [N]
    real(8),dimension(0:2)::Force !機体にはたらく外力 [X,Y,Z] [N]
    real(8),dimension(0:2)::Moment !機体にはたらくモーメント [L,M,N] [N*m]
    real(8)::epsilon !主翼位置での吹きおろし角 [deg]
    !機体データ
    integer::span_div,chord_div !スパン方向分割数
    real(8)::Mass !機体重量 [kg]
    real(8)::Ixx,Iyy,Izz,Ixz !慣性モーメント [kg*m^4]
    real(8)::CD_fuselage,S_fuselage,chord_fuselage !カウルの抗力係数,面積，コード長
    !その他
    integer::pilot !操縦方法(0:重心移動，1:エレベーター)
    integer::Tail_configulation !尾翼配置(0:通常，1:T字尾翼)
    real(8)::Drag_fuselage
    real(8)::CL,CD,Cm,CL_old
    real(8)::rho,mu,Re
    real(8)::Velocity
    real(8)::X,Y,Z,L,M,N
    real(8)::phi,theta,psi
    real(8)::de,dh,dr
    real(8)::p,q,r
    real(8)::u,v,w,ug,vg,wg
    real(8)::alpha,beta,gamma
    real(8)::chord_mac,at !MAC,水平尾翼揚力傾斜(オールフライング，対称翼を仮定) [N/deg]
    real(8)::lt,it,Lift_tail,M0_wing,Lift_wing,hac,hcg
    real(8)::pi,deg,rad !円周率，degへの変換，radへの変換
    real(8)::tmp
    !カウンター
    integer::i,j,k
    !計算時間計測
    real(8)::CPU_time_s
    integer::initial_time,finish_time,time_rate,CPU_time_min !initial time, finish time, time rate

    pi=4.0D0*atan(1.0D0) !円周率
    deg=(180.0D0/pi) !degへの変換
    rad=(pi/180.0D0) !radへの変換

    Altitude_angle=0.0D0 !姿勢角 [Φ,θ,Ψ] [deg]
    Ground_speed=0.0D0 !機体軸における対地速度ベクトル [u,v,w] [m/s]
    Angular_velocity=0.0D0 !機体軸における角速度ベクトル [p,q,r] [deg/s]
    Input=0.0D0 !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg] [δe,δh,δr]
    gust=0.0D0 !機体軸における突風ベクトル [ug,vg,wg] [m/s]
    Acceleration=0.0D0 !機体軸における対地加速度ベクトル [du/dt,dv/dt,dw/dt] [m/s]
    Angular_acceleration=0.0D0 !機体軸における角加速度ベクトル [dp/dt,dq/dt,dr/dt] [deg/s]

    !機体データ読み込み
    open(10,file="Plane_data.txt") !ファイルを開く
        read(10,*) Mass !機体重量 [kg]
        read(10,*) Ixx,Iyy,Izz,Ixz !慣性モーメント [kg*m^4]
        read(10,*) CD_fuselage,S_fuselage,chord_fuselage !カウルの抗力係数,面積，コード長
    close(10) !ファイルを閉じる

    !わかりやすいように変数変換
    u=Ground_speed(0);          w=Ground_speed(2);              v=Ground_speed(1)           !u,v,w [m/s]
    ug=gust(0);                 wg=gust(2);                     vg=gust(1)                  !u,v,w [m/s]
    phi=Altitude_angle(0);      theta=Altitude_angle(1);            psi=Altitude_angle(2)   !Φ，θ，ψ [rad]を格納
    p=Angular_velocity(0);      q=Angular_velocity(1);          r=Angular_velocity(2)       !p,q,r [deg]       
    de=Input(0);                dh=Input(1);                    dr=Input(2)                 !操舵を格納

    !速度V [m/s]，全機迎角α [deg]，横滑り角β [deg]に初期値を入れる
    beta=0.0D0

    write(*,*) "How pilot?"
    write(*,*) "0: dh"
    write(*,*) "1: de"
    read(*,*) pilot
    write(*,*) "Tail Configulation?"
    write(*,*) "0: Conventional"
    write(*,*) "1: T-tail"
    read(*,*) Tail_configulation
    write(*,*) "hE?"
    read(*,*) hE
    write(*,*)
    write(*,'(A12,A12,A12,A12,A12,A12,A12,A12,A12,A12)') &
        "alpha [deg]","V [m/s]","CL [-]","CD [-]","Cm [-]","gamma [deg]","u [m/s]","w [m/s]","de [deg]","dh [deg]"
    
    call system_clock(initial_time) !開始時間の読み込み

    CL_old=0.0D0 !初期値0
    k=0
    do i=0,int((alpha_max-alpha_min)/da)
            
        alpha=alpha_min+da*real(i)
        Velocity=Velocity_max
        !write(*,'(F12.3)') alpha
        
        do j=0,10000

            Lift=0.0D0 !揚力 [N]
            Drag=0.0D0 !抗力 [N]
            Force=0.0D0 !機体にはたらく外力 [X,Y,Z] [N]
            Moment=0.0D0 !機体にはたらくモーメント [L,M,N] [N*m]

            !主翼計算    
            open(10,file="Wing_data.txt") !ファイルを開く
            read(10,*) span_div !スパン方向分割数，コード方向分割数
            read(10,*) rho,mu !空気密度 [kg/m^3]，動粘性係数 [m^2/s]
            read(10,*) hac,tmp !空力中心位置
            read(10,*) tmp !翼素間隔 [m]
            read(10,*) chord_mac,hcg !平均空力翼弦長 [m]，重心位置 [-]
            close(10) !ファイルを閉じる
            !call Wing_calculation(Lift,Drag,Force,Moment,Velocity,alpha,beta,p,r,hE,dh,epsilon,span_div)
            call Wing_calculation_linear(Lift,Drag,Force,Moment,Velocity,alpha,beta,p,r,hE,dh,epsilon,span_div)
            if(isnan(Lift)) write(*,'(A24)') "Wing"
            Lift_wing=Lift
            M0_wing=Moment(1)-Lift_wing*chord_mac*(hcg-dh-hac)
            if(Tail_configulation==1) epsilon=0.0D0

            !水平尾翼計算
            open(10,file="Tail_data.txt") !ファイルを開く
            read(10,*) span_div,chord_div !スパン方向分割数，コード方向分割数
            read(10,*) tmp,tmp !空気密度 [kg/m^3]，動粘性係数 [m^2/s]
            read(10,*) tmp,tmp !空力中心位置
            read(10,*) tmp !翼素間隔 [m]
            read(10,*) tmp,lt,it !平均空力翼弦長 [m]，重心-水平尾翼中心間距離 [m]，水平尾翼取り付け角 [deg]
            close(10) !ファイルを閉じる
            call Tail_calculation(Lift,Drag,Force,Moment,Velocity,alpha,q,de,dh,epsilon,span_div,chord_div)
            if(isnan(Lift)) write(*,'(A24)') "Tail"
            Lift_tail=Lift-Lift_wing !尾翼揚力を計算 [N]
            at=Lift_tail/(it-epsilon) !尾翼揚力傾斜を計算 [N/deg]

            !垂直尾翼計算
            open(10,file="Fin_data.txt") !ファイルを開く
                read(10,*) span_div,chord_div !スパン方向分割数読み込み
            close(10) !ファイルを閉じる
            call Fin_calculation(Lift,Drag,Force,Moment,Velocity,beta,p,r,dr,dh,span_div,chord_div)

            !カウル抗力計算
            Drag_fuselage=0.5D0*rho*Velocity**2*CD_fuselage*S_fuselage
            Drag=Drag+Drag_fuselage

            !CL,CD,Cm,Gammaを計算
            CL=Lift/(0.5D0*rho*Velocity**2*S)
            CD=Drag/(0.5D0*rho*Velocity**2*S)
            Cm=Moment(1)/(0.5D0*rho*Velocity**2*S*chord_mac)
            gamma=-deg*atan(Drag/Lift)
            u=Velocity*cos(rad*gamma)
            w=Velocity*sin(rad*gamma)

            if(pilot==0) then !M=0になるように重心位置を変える
                if(Cm>0.0D0) dh=dh+0.0005
                if(Cm<0.0D0) dh=dh-0.0005
            end if
            if(pilot==1) then !M=0になるようにエレベータを操作する
                if(Cm>0.0D0) de=de+0.001
                if(Cm<0.0D0) de=de-0.001
            end if

            !if(abs((Lift-Mass*g)/Mass*g)<error .and. abs(Cm)<error) Then
            if(abs(Cm)<error) Then
                write(*,'(F12.6,F12.6,F12.6,F12.6,F12.6,F12.6,F12.6,F12.6,F12.6,F12.6)') alpha,Velocity,CL,CD,Cm,gamma,u,w,de,dh
                exit
            elseif(Velocity>Velocity_max) then
                exit
            else
                Velocity=sqrt((2.0D0*Mass*g)/(rho*S*CL))
            end if
        end do 

        if(CL<CL_old) then !失速したら終了
            if(k>0) exit
            k=k+1
        else
            CL_old=CL
        end if

    end do

    call system_clock(finish_time,time_rate) !終了時間の読み込み
    CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
    if(CPU_time_s<60) then
        write(*,'(A24,f12.3,A12)') "Computation Time",CPU_time_s,"[s]"
    else 
        CPU_time_min=int(CPU_time_s/60.0D0)
        CPU_time_s=CPU_time_s-60.D0*CPU_time_min
        write(*,'(A24,I12,A12,F12.3,A12)') "Computation Time",CPU_time_min,"[min]",CPU_time_s,"[s]"
    end if

end subroutine

subroutine Test_NLFS()
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan !プログラム中で使用する関数の名前を宣言
    real(8)::Distance !飛行距離 [m]
    integer,parameter::genom_length = 21 !遺伝子情報の長さ
    real(8),dimension(0:genom_length-1)::genom_list !遺伝子情報
    integer,parameter::FLightLog_data = 47 !飛行解析のデータ数
    real(8),dimension(:,:),allocatable::FlightLog !飛行解析のデータ
    real(8)::time_max !最大飛行時間 [s]
    real(8)::time_step !タイムステップ [s]
    integer::iteration_max !最大反復回数 [s]
    integer::i_max_copy !最大反復回数のコピー [s]
    integer::pilot_method !操縦方式．0：エレベータ，1：重心移動
    !カウンター
    integer::i,j
    !計算時間計測
    real(8)::CPU_time_s
    integer::initial_time,finish_time,time_rate,CPU_time_min,CPU_time_h !initial time, finish time, time rate
    
    call system_clock(initial_time) !開始時間の読み込み

    !遺伝子情報の読み込み
    open(10,file="LIST.txt") !ファイルを開く
    read(10,*) genom_list(0),genom_list(1),genom_list(2)
    do i=3,8
        read(10,*) genom_list(i) !1行ずつ値を読み込む
    end do
    do i=0,3
        read(10,*) genom_list(9+i*3),genom_list(10+i*3),genom_list(11+i*3)
    end do
    close(10) !ファイルを閉じる

    write(*,*) "-----------------------------------------------------------------------------------"
    write(*,'(A15)') "GENOM"
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "phi","=",genom_list(0),"theta","=",genom_list(1),"psi","=",genom_list(2)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "time_pullup","=",genom_list(3),"dive_angle","=",genom_list(6),"cluse_angle","=",genom_list(7)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "time_bank","=",genom_list(4),"time_bank_end","=",genom_list(5),"bank_angle","=",genom_list(8)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "PID_lon1[0]","=",genom_list(9),"PID_lon2[0]","=",genom_list(12)&
    ,"PID_lat1[0]","=",genom_list(15),"PID_lat2[0]","=",genom_list(18)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "PID_lon1[1]","=",genom_list(10),"PID_lon2[1]","=",genom_list(13)&
    ,"PID_lat1[1]","=",genom_list(16),"PID_lat2[1]","=",genom_list(19)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "PID_lon1[1]","=",genom_list(11),"PID_lon2[2]","=",genom_list(14)&
    ,"PID_lat1[2]","=",genom_list(17),"PID_lat2[2]","=",genom_list(20)
    write(*,*)

    !飛行時間，タイムステップの読み込み
    open(10,file="NLFS.txt") !ファイルを開く
    read(10,*) pilot_method !操縦方式
    read(10,*) time_max !最大飛行時間
    read(10,*) time_step !タイムステップ
    close(10) !ファイルを閉じる
    !call write1D(genom_length,genom_list)

    iteration_max=int(time_max/time_step)+1 !最大反復回数の読み込み
    allocate(FlightLog(0:iteration_max,0:FLightLog_data))
    write(*,*) iteration_max 

    i_max_copy=iteration_max
    call NLFS(Distance,genom_length,genom_list,i_max_copy,pilot_method,FlightLog_data,FlightLog)

    !call write2D(iteration_max+1,FlightLog_data+1,FlightLog)
    
    call Output_FlightLog(i_max_copy,iteration_max,FlightLog_data,FlightLog,genom_length,genom_list) !飛行解析データをcsvファイルで書き出し

    write(*,*) 
    write(*,'("Flight Distance =",F8.3,"[m]")') Distance
    write(*,*) 
    
    write(*,*) "------------------------------------------------------------"
    call system_clock(finish_time,time_rate) !終了時間の読み込み
    CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
    if(CPU_time_s>60) then
        CPU_time_min=int(CPU_time_s/60.0D0)
        CPU_time_s=CPU_time_s-60.D0*CPU_time_min
        if(CPU_time_min>60) then
            CPU_time_h=int(CPU_time_min/60)
            CPU_time_min=CPU_time_min-60*CPU_time_h
        end if
    end if
    write(*,'(A24,I12,A12,I12,A12,F12.3,A12)') "Computation Time",CPU_time_h,"[h]",CPU_time_min,"[min]",CPU_time_s,"[s]"

end subroutine

subroutine NLFS(Distance,genom_length,genom_list,iteration_max,pilot_method,FlightLog_data,FlightLog)
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan,sqrt,sin,cos,tan,mod !プログラム中で使用する関数の名前を宣言
    !パラメーター
    real(8),parameter::relaxation_coef =1.00D0 !操舵の緩和係数．舵角の急激な変化を抑える
    real(8),parameter::weighting_factor1 =0.000D0 !飛距離と偏差の総和の重みづけ．大きければ大きいほど偏差の総和の少なさが優先される．2.5
    real(8),parameter::weighting_factor2 =1.00D0 !縦と横・方向の偏差の総和の重みづけ．大きければ大きいほど横・方向の少なさが優先される．10
    real(8),parameter::weighting_factor3 =0.000D0 !飛距離に対する入力の時間変化率の総和の重みづけ．大きければ大きいほど入力が滑らかになることが優先される．
    real(8),parameter::gravity=9.795 !重力加速度 [m/s^2]
    !引数
    real(8),intent(OUT)::Distance !飛行距離
    integer,intent(IN)::genom_length !遺伝子情報の長さ
    real(8),dimension(0:genom_length-1),intent(IN)::genom_list !遺伝子情報
    integer,intent(INOUT)::iteration_max !最大反復回数
    integer,intent(INOUT)::pilot_method !操縦方式．0：エレベータ，1：重心移動
    integer,intent(IN)::FlightLog_data !FLightLogのデータ数
    real(8),dimension(0:iteration_max,0:FlightLog_data),intent(OUT)::FlightLog !飛行解析のデータ
    !非線形解析用
    real(8)::Lift !揚力 [N]
    real(8)::Drag !抗力 [N]
    real(8),dimension(0:2)::Force !機体にはたらく外力 [X,Y,Z] [N]
    real(8),dimension(0:2)::Moment !機体にはたらくモーメント [L,M,N] [N*m]
    real(8),dimension(0:2)::Position !位置ベクトル [xE,yE,hE] [m]
    real(8),dimension(0:2)::Velocity !速度ベクトル [dxE/dt,dyE/dt,dhE/dt] [m/s]
    real(8),dimension(0:2)::Velocity_old !速度ベクトル [dxE/dt,dyE/dt,dhE/dt] [m/s]
    real(8),dimension(0:2)::Altitude_angle !姿勢角 [Φ,θ,Ψ] [deg]
    real(8),dimension(0:2)::Altitude_angular_velocity !姿勢角速度 [dΦ/dt,dθ/dt,dΨ/dt] [deg]
    real(8),dimension(0:2)::Alt_ang_vel_old !姿勢角速度 [dΦ/dt,dθ/dt,dΨ/dt] [deg]
    real(8),dimension(0:2)::Ground_speed !機体軸における対地速度ベクトル [u,v,w] [m/s]
    real(8),dimension(0:2)::Acceleration !機体軸における対地加速度ベクトル [du/dt,dv/dt,dw/dt] [m/s]
    real(8),dimension(0:2)::Acceleration_old !機体軸における対地加速度ベクトル [du/dt,dv/dt,dw/dt] [m/s]
    real(8),dimension(0:2)::Angular_velocity !機体軸における角速度ベクトル [p,q,r] [deg/s]
    real(8),dimension(0:2)::Angular_acceleration !機体軸における角加速度ベクトル [dp/dt,dq/dt,dr/dt] [deg/s]
    real(8),dimension(0:2)::Angular_acceleration_old !機体軸における角加速度ベクトル [dp/dt,dq/dt,dr/dt] [deg/s]
    real(8),dimension(0:2)::Air_angle !迎角，横滑り角，経路角 [α,β,γ] [deg]
    real(8),dimension(0:2)::Input !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg] [δe,δh,δr]
    real(8),dimension(0:2)::Input_old !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg] [δe,δh,δr]
    real(8),dimension(0:2)::Wind !対地風ベクトル [ugE,vgE,wgE] [m/s] @プラホ高度における
    real(8),dimension(0:2)::gust !機体軸における突風ベクトル [ug,vg,wg] [m/s]
    real(8),dimension(0:1)::Speed !対地速度，対気速度 [VE,V] [m/s]
    real(8),dimension(0:2)::Trim_speed !トリム速度ベクトル [u0,v0,w0] [m/s]
    real(8),dimension(0:2)::Trim_angle !トリム角度 [Φ0,θ0,Ψ0] [m/s]
    real(8),dimension(0:2)::Target !目標値．ダイブ角 [deg]，定常滑空角 [deg]，バンク角 [deg]
    real(8),dimension(0:2)::PID_para_lon1,PID_para_lon2 !縦の制御パラメーター(ダイブ),(定常滑空) [-]
    real(8),dimension(0:2)::PID_para_lat1,PID_para_lat2 !横・方向の制御パラメーター(定常滑空),(バンク) [m/s]
    real(8),dimension(0:2)::PID_lon1,PID_lon2 !縦の制御パラメーター(ダイブ)，(定常滑空)
    real(8),dimension(0:2)::PID_lat1,PID_lat2 !横・方向の制御パラメーター(定常滑空),(バンク)
    real(8),dimension(0:iteration_max)::epsilon_history !主翼位置での吹きおろしε [deg] の時間履歴
    real(8)::Trim_velocity !トリム速度 V [m/s]
    real(8)::epsilon !主翼位置での吹きおろしε [deg]
    real(8)::lt !重心-水平尾翼空力中心間距離 [m]
    real(8)::Mass,Weight !機体質量 [kg]，機体重量 [N]
    real(8)::nx,ny,nz !荷重倍数
    !着水判定用
    real(8)::alpha_stall !失速角 [deg]
    real(8)::hE_water !着水高度
    real(8),parameter::limit_load_factor =3.0D0 !制限荷重倍数 [G]
    !操舵の制限
    real(8)::de_max,dh_max,dr_max !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg]
    !評価関数
    real(8)::sum_error,sum_D_input !偏差の総和，入力の時間変化率の総和
    !時間
    real(8)::time !時間 [s]
    real(8)::time_max !最大飛行時間 [s]
    real(8)::time_step !タイムステップ [s]
    real(8)::time_pullup !引き起こし終了時刻 [s]
    real(8)::time_bank !バンク開始時刻 [s]
    real(8)::time_bank_end !バンク終了時刻 [s]
    !その他
    real(8)::phi,theta,psi
    real(8)::p,q,r
    real(8)::u,v,w,ug,vg,wg
    real(8)::xE,yE,hE,hE_0
    real(8)::alpha,beta,gamma
    real(8)::tmp,de,dh,dr
    real(8),dimension(:,:),allocatable::matrix,matrix1,matrix2,matrix3
    real(8)::pi,deg,rad !円周率，degへの変換，radへの変換
    !カウンター
    integer::i,j,k,l,m,n
    integer::iteration
    !計算時間計測
    real(8)::CPU_time_s
    integer::initial_time,finish_time,time_rate,CPU_time_min !initial time, finish time, time rate

    pi=4.0D0*atan(1.0D0) !円周率
    deg=(180.0D0/pi) !degへの変換
    rad=(pi/180.0D0) !radへの変換
    
    call system_clock(initial_time) !開始時間の読み込み

    !配列の初期化
    Lift=0.0D0;         Drag=0.0D0
    Force=0.0D0;        Moment=0.0D0
    Position=0.0D0;     Velocity=0.0D0; Velocity_old=0.0D0 !速度ベクトル [m]
    Altitude_angle=0.0D0;               Altitude_angular_velocity=0.0D0 !姿勢角 [deg]
    Alt_ang_vel_old=0.0D0 !姿勢角 [deg]
    Ground_speed=0.0D0 !機体軸における対地速度ベクトル [m/s]
    Acceleration=0.0D0;                 Acceleration_old=0.0D0
    Angular_velocity=0.0D0;             Angular_acceleration=0.0D0
    Angular_acceleration_old=0.0D0
    Air_angle=0.0D0;    Input=0.0D0;    Input_old=0.0D0 !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg]
    Wind=0.0D0;         gust=0.0D0;     speed=0.0D0
    PID_lon1=0.0D0;     PID_lon2=0.0D0; PID_lat1=0.0D0;     PID_lat2=0.0D0
    epsilon_history=0.0D0;              epsilon=0.0D0
    FlightLog=0.0D0;    sum_error=0.0D0;sum_D_input=0.0D0;  Distance=0.0D0
    CPU_time_min=0.0D0; time=0.0D0

    !---------------------------------------------------!
    !--------------------値の読み込み--------------------!
    !---------------------------------------------------!

    do i=0,2
        Altitude_angle(i)=genom_list(i)
    end do
    time_pullup=genom_list(3) !引き起こし時間 [s]
    time_bank=genom_list(4) !バンク開始時間 [s]
    time_bank_end=genom_list(5) !バンク開始時間 [s]
    do i=0,2
        Target(i)=genom_list(6+i)
        PID_para_lon1(i)=genom_list(9+i)
        PID_para_lon2(i)=genom_list(12+i)
        PID_para_lat1(i)=genom_list(15+i)
        PID_para_lat2(i)=genom_list(18+i)
    end do
    !write(*,'(A24,I12)') "iteration_max",iteration_max
    !write(*,*) time_pullup,time_bank,time_bank_end
    !call write1D(3,PID_para_lon1)
    
    write(*,*) "------------------------------------------------------------"
    !write(*,*)
    write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') &
    "phi",genom_list(0),"theta",genom_list(1),"psi",genom_list(2)
    write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') &
    "pullup",genom_list(3),"dive",genom_list(6),"cluse",genom_list(7)
    write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') &
    "bank",genom_list(4),"to",genom_list(5),"angle",genom_list(8)
    write(*,*)

    open(10,file="NLFS.txt") !ファイルを開く
        read(10,*) pilot_method !操縦方式
        read(10,*) time_max !最大飛行時間 [s]
        read(10,*) time_step !タイムステップ [s]
        read(10,*) alpha_stall !失速角 [s]
        read(10,*) hE_water !着水高度 [s]
        !write(*,*) time_max,time_step,alpha_stall,hE_water
        read(10,*) Wind(0:2) !対地風ベクトル [m/s]
        read(10,*) Trim_velocity !トリム速度V [m/s]
        !トリム速度
        read(10,*) Trim_speed(0:2)
        !トリム姿勢角
        read(10,*) Trim_angle(0:2)
        !初期条件
        read(10,*) Position(0:2) !位置
        read(10,*) speed(0) !プラホからの飛び出し速度 [m/s]
        read(10,*) Air_angle(2) !経路角 [deg] プラホなら -4 [deg]
        !舵角の制限
        read(10,*) de_max
        read(10,*) dh_max
        read(10,*) dr_max
    close(10) !ファイルを閉じる
    hE_0=Position(2) !プラホ高度 [m]

    !重心-水平尾翼空力中心間距離の読み込み
    open(10,file="Tail_data.txt") !ファイルを開く
        do i=0,3
            read(10,*) tmp
        end do
        read(10,*) tmp,lt,tmp ! ，重心-水平尾翼中心間距離 [m] 
    close(10) !ファイルを閉じる
    !write(*,*) lt

    !機体データ読み込み
    open(10,file="Plane_data.txt") !ファイルを開く
        read(10,*) Mass !機体重量 [kg]
    close(10) !ファイルを閉じる
    Weight=Mass*gravity !機体重量を計算 [N]

    !---------------------------------------------------!
    !--------------------初期値の計算--------------------!
    !---------------------------------------------------!

    !Φ，θ，ψ [rad]を格納
    phi=rad*Altitude_angle(0);    theta=rad*altitude_angle(1);    psi=rad*altitude_angle(2)

    !ug,vg,wg [m/s]を計算
    allocate(matrix1(0:2,0:0),matrix2(0:2,0:0))
    do i=0,2
        matrix1(i,0)=wind(i) !対地風ベクトルをmatrix1に格納
    end do
    call Transform_coordinate(phi,theta,psi,matrix1,matrix2) !絶対座標系から機体軸座標に変換する
    !call Transform_coordinate(phi,Trim_angle(1),psi,matrix1,matrix2) !絶対座標系から機体軸座標に変換する
    do i=0,2
        gust(i)=matrix2(i,0)
    end do
    deallocate(matrix1,matrix2)
    ug=gust(0);     vg=gust(1);     wg=gust(2)

    !u,v,w [m/s] ug,vg,wg [m/s]を計算
    Ground_speed(0)=speed(0)*cos(theta-Air_angle(2)*(pi/180))
    Ground_speed(1)=0.0D0
    Ground_speed(2)=speed(0)*sin(theta-Air_angle(2)*(pi/180))
    u=ground_speed(0);  v=ground_speed(1);  w=ground_speed(2);  

    !対地速度VEの計算
    speed(1)=sqrt((u+ug)*(u+ug)+(v+vg)*(v+vg)+(w+wg)*(w+wg))

    !迎角，横滑り角，経路角の初期値の計算
    Air_angle(0)=(180.0D0/pi)*atan(w/u)
    Air_angle(1)=(180.0D0/pi)*atan(v/speed(0))
    Air_angle(2)=Altitude_angle(1)-Air_angle(0)

    !---------------------------------------------------!
    !--------------------反復計算開始--------------------!
    !---------------------------------------------------!

    do iteration=0,iteration_max

        u=Ground_speed(0);          w=Ground_speed(2);              v=Ground_speed(1)           !u,v,w [m/s]
        phi=rad*Altitude_angle(0);  theta=rad*altitude_angle(1);    psi=rad*altitude_angle(2)   !Φ，θ，ψ [rad]を格納
        p=Angular_velocity(0);      q=Angular_velocity(1);          r=Angular_velocity(2)       !p,q,r [deg]
        xE=Position(0);             yE=Position(1);                 hE=Position(2)              !位置xE,yE,hE [m]を格納       
        alpha=Air_angle(0);         beta=Air_angle(1);              gamma=air_angle(2)          !迎角α,横滑り角β [deg] を格納       
        de=Input(0);                dh=Input(1);                    dr=Input(2)                 !操舵を格納

        call system_clock(finish_time,time_rate) !終了時間の読み込み
        CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
        if(CPU_time_s>60) then
            CPU_time_min=int(CPU_time_s/60.0D0)
            CPU_time_s=CPU_time_s-60.D0*CPU_time_min
        end if
        if(mod(iteration-1,1)==0 .or. iteration==0) then
        !if(iteration<0) then
        !    write(*,*)
        !    write(*,'(A12,F12.3)') "time",time
        !    write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "xE",xE,"yE",yE,"hE",hE
        !    write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "u",u,"v",v,"w",w
        !    write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "phi",deg*phi,"theta",deg*theta,"psi",deg*psi
        !    write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "p",p,"q",q,"r",r
        !    write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "alpha",alpha,"beta",beta,"gamma",gamma
        !    write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "de",de,"dh",dh,"dr",dr
        !    write(*,'(A24,I12,A12,F12.3,A12)') "Computation Time",CPU_time_min,"[min]",CPU_time_s,"[s]"
        !    write(*,'(A12,F12.3,A12,F12.3,A12,F12.3,A12,F12.3)') "time",time,"xE",xE,"yE",yE,"hE",hE
        end if
        

        !-------------------------------------------------!
        !--------------------時間の更新--------------------!
        !-------------------------------------------------!

        time=time_step*iteration
        !write(*,'(A12,F12.3)') "time",time
        
        !--------------------------------------------------------------!
        !--------------------FlightLogにデータを入力--------------------!
        !--------------------------------------------------------------!

        FlightLog(iteration,0)=time
        do i=0,2
            FlightLog(iteration,1+i)=Position(i) !位置ベクトル [m]
            FlightLog(iteration,4+i)=Velocity(i) !速度ベクトル [m/s]
            FlightLog(iteration,7+i)=Altitude_angle(i) !姿勢角 [deg]
            FlightLog(iteration,10+i)=Altitude_angular_velocity(i) !姿勢角速度 [deg]
            FlightLog(iteration,13+i)=Ground_speed(i) !機体軸における対地速度ベクトル [m/s]
            FlightLog(iteration,16+i)=Angular_velocity(i) !機体軸における角速度ベクトル [deg/s]
            FlightLog(iteration,19+i)=Air_angle(i) !迎角，横滑り角，経路角 [deg]
            FlightLog(iteration,22+i)=Input(i) !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg]
            FlightLog(iteration,25+i)=Wind(i) !対地風ベクトル [m/s]
            FlightLog(iteration,28+i)=gust(i) !突風ベクトル [m/s]
        end do
        do i=0,1
            FlightLog(iteration,31+i)=speed(i) !対地速度，対気速度 [m/s]
        end do
        FlightLog(iteration,33)=Lift !揚力 [N]
        FlightLog(iteration,34)=Drag !抗力 [N]
        do i=0,2
            FlightLog(iteration,35+i)=Force(i) !機体にはたらく外力 [X,Y,Z] [N]
            FlightLog(iteration,38+i)=Moment(i) !機体にはたらくモーメント [L,M,N] [N*m]
        end do
        FlightLog(iteration,41)=nx !荷重倍数 [-]
        FlightLog(iteration,42)=ny !荷重倍数 [-]
        FlightLog(iteration,43)=nz !荷重倍数 [-]
        FlightLog(iteration,44)=Lift/Drag !揚抗比
        FlightLog(iteration,45)=sqrt(ug*ug+vg*vg+wg*wg) !現在の高度での風速
        if(time<time_pullup) then !ダイブ中
            FlightLog(iteration,46)=Target(0) !目標経路角
        else !定常滑空
            FlightLog(iteration,46)=Target(1) !目標経路角
        end if
        if(time<time_bank .or. time>time_bank_end) then !定常滑空(目標バンク角0 [deg])
            FlightLog(iteration,47)=0 !目標バンク角
        elseif(time>time_bank .and. time<time_bank_end) then !バンク中
            FlightLog(iteration,47)=Target(2) !目標バンク角
        end if

        !-----------------------------------------------------------------!
        !--------------------対地風ベクトルを機体軸に変換--------------------!
        !-----------------------------------------------------------------!
        
        !対地風ベクトルをmatrix1に格納
        allocate(matrix1(0:2,0:0),matrix2(0:2,0:0))
        do i=0,2
            matrix1(i,0)=wind(i)
        end do
        call Transform_coordinate(phi,theta,psi,matrix1,matrix2) !絶対座標系から機体軸座標に変換する
        do i=0,2
            gust(i)=matrix2(i,0)
        end do
        deallocate(matrix1,matrix2)
        !1/6乗則によって地面付近の風を弱くする．
        gust=gust*((hE/hE_0)**(1.0D0/6.0D0))

        !-----------------------------------------------!
        !--------------------PID制御--------------------!
        !-----------------------------------------------!

        !P,I,Dを計算
        !縦
        if(time<time_pullup) then !ダイブ中
            PID_lon1(2)=((Air_angle(2)-Target(0))-PID_lon1(0))/time_step ![deg/s]
            PID_lon1(0)=(Air_angle(2)-Target(0)) ![deg]
            PID_lon1(1)=PID_lon1(1)+PID_lon1(0)*time_step ![deg*s]
        else !定常滑空
            PID_lon2(2)=((Air_angle(2)-Target(1))-PID_lon2(0))/time_step ![deg/s]
            PID_lon2(0)=(Air_angle(2)-Target(1)) ![deg]
            PID_lon2(1)=PID_lon2(1)+PID_lon2(0)*time_step ![deg*s]
        end if
        !横・方向
        if(time<time_bank .or. time>time_bank_end) then !定常滑空(目標バンク角0 [deg])
            PID_lat1(2)=((Altitude_angle(0)-0.0D0)-PID_lat1(0))/time_step ![deg/s]
            PID_lat1(0)=(Altitude_angle(0)-0.0D0) ![deg]
            PID_lat1(1)=PID_lat1(1)+PID_lat1(0)*time_step ![deg*s]
        elseif(time>time_bank .and. time<time_bank_end) then !バンク中
            PID_lat2(2)=((Altitude_angle(0)-Target(2))-PID_lat2(0))/time_step ![deg/s]
            PID_lat2(0)=(Altitude_angle(0)-Target(2)) ![deg]
            PID_lat2(1)=PID_lat2(1)+PID_lat2(0)*time_step ![deg*s]
        end if

        !操舵量を計算．エレベータ舵角も重心移動距離も負のピッチングモーメントを発生させる向きが正．
        !縦
        Input=0.0D0
        if(time<time_pullup) then !ダイブ中
            do i=0,2
                Input(pilot_method)=Input(pilot_method)+PID_para_lon1(i)*PID_lon1(i)
            end do
        else !定常滑空
            do i=0,2
                Input(pilot_method)=Input(pilot_method)+PID_para_lon2(i)*PID_lon2(i)
            end do
        end if
        IF(pilot_method==0) then !エレベータで操縦するなら
            Input(0)=Input(0)*(de_max/dh_max) !PIDの値を調整
        end if
        !横・方向
        if(time<time_bank .or. time>time_bank_end) then !定常滑空(目標バンク角0 [deg])
            do i=0,2
                Input(2)=Input(2)+PID_para_lat1(i)*PID_lat1(i)*(dr_max/dh_max)
            end do
        elseif(time>time_bank .and. time<time_bank_end) then !バンク中
            do i=0,2
                Input(2)=Input(2)+PID_para_lat2(i)*PID_lat2(i)*(dr_max/dh_max)
            end do
        end if
        !緩和係数をかける
        do i=0,2
            Input(i)=Input_old(i)+relaxation_coef*(Input(i)-input_old(i))
        end do 
        !舵角を制限
        if(abs(Input(0))>de_max) then
            Input(0)=(Input(0)/abs(Input(0)))*de_max
        end if
        if(abs(Input(1))>dh_max) then
            Input(1)=(Input(1)/abs(Input(1)))*dh_max
        end if
        if(abs(input(2))>dh_max) then
            Input(2)=(Input(2)/abs(Input(2)))*dr_max
        end if
        !Input_oldの更新はこのループの最後で

        !------------------------------------------------!
        !--------------------全機計算--------------------!
        !------------------------------------------------!

        !尾翼位置での吹きおろしを計算
        if(iteration-int(lt/(speed(1)*time_step))>0) then
            epsilon=epsilon_history(iteration-int(lt/(speed(1)*time_step)))
        else
            epsilon=0.000D0
        end if

        !write(*,'(A12,F12.3)') "epsilon",epsilon
        call Plane_calculation &
            (Altitude_angle,Ground_speed,Angular_velocity,gust,Input &
            ,hE,Lift,Drag,Force,Moment,epsilon,Acceleration,Angular_acceleration)

        epsilon_history(iteration)=epsilon

        !---------------------------------------------------!
        !--------------------加速度を積分--------------------!
        !---------------------------------------------------!

        !Adams-Bashforth法で積分
        Ground_speed=Ground_speed+((3.0D0*Acceleration-Acceleration_old)/2.0D0)*time_step
        Angular_velocity=Angular_velocity+((3.0D0*Angular_acceleration-Angular_acceleration_old)/2.0D0)*time_step

        !oldを更新
        Acceleration_old=Acceleration
        Angular_acceleration_old=Angular_acceleration

        !---------------------------------------------------!
        !--------------------角速度を計算--------------------!
        !---------------------------------------------------!
        
        p=Angular_velocity(0)
        q=Angular_velocity(1)
        r=Angular_velocity(2)
        
        !機体軸p,q,rから慣性系dΦ/dt,dθ/dt,dΨ/dtへ
        Altitude_angular_velocity(0)=p+(r*cos(phi)+q*sin(phi))*tan(theta)
        Altitude_angular_velocity(1)=q*cos(phi)-r*sin(phi)
        Altitude_angular_velocity(2)=(r*cos(phi)+q*sin(phi))/cos(theta)

        !-------------------------------------------------!
        !--------------------速度を計算--------------------!
        !-------------------------------------------------!
        
        !速度を計算
        ug=gust(0)
        vg=gust(1)
        wg=gust(2)
        u=Ground_speed(0)
        w=Ground_speed(2)
        v=Ground_speed(1)
        speed(0)=sqrt(u*u+v*v+w*w)
        speed(1)=sqrt((u+ug)**2+(v+vg)**2+(w+wg)**2)

        !機体軸u,v,wから慣性系dxE/dt,dyE/dt,dhE/dtへ
        allocate(matrix1(0:2,0:0),matrix2(0:2,0:0))
        do i=0,2
            matrix1(i,0)=Ground_speed(i)
        end do
        call Inverse_Transform_coordinate(phi,theta,psi,matrix1,matrix2) !機体座標系から絶対軸座標に変換する
        do i=0,2
            Velocity(i)=matrix2(i,0)
        end do
        Velocity(2)=-Velocity(2) !dhE/dtは上向きが正
        deallocate(matrix1,matrix2)

        !-----------------------------------------------------------------!
        !--------------------位置，姿勢角,荷重倍数を計算--------------------!
        !-----------------------------------------------------------------!

        !Adams-Bashforth法で積分
        Altitude_angle=Altitude_angle+((3.0D0*Altitude_angular_velocity-Alt_ang_vel_old)/2.0D0)*time_step
        Position=Position+((3.0D0*Velocity-Velocity_old)/2.0D0)*time_step
        
        !dx/dt_oldを更新
        Alt_ang_vel_old=Altitude_angular_velocity
        Velocity_old=Velocity
        
        !迎角，横滑り角，経路角の初期値の計算
        Air_angle(0)=(180.0D0/pi)*atan((w+wg)/(u+ug))
        Air_angle(1)=(180.0D0/pi)*atan((v+vg)/speed(0))
        Air_angle(2)=Altitude_angle(1)-Air_angle(0)

        !荷重倍数の計算
        nx=Force(0)/Weight
        ny=Force(1)/Weight
        nz=-Force(2)/Weight

        !-----------------------------------------------!
        !--------------------着水判定--------------------!
        !-----------------------------------------------!

        tmp=0.0D0
        if(iteration==iteration_max) then
            write(*,'(A24)') "Over Flight time"
            tmp=tmp+1.0D0
        elseif(alpha>alpha_stall) then
            write(*,'(A24,F12.5,F12.5)') "Stalled",alpha,alpha_stall
            tmp=tmp+1.0D0
        elseif(position(2)<hE_water) then
            write(*,'(A24,F12.5,F12.5)') "Landing",Position(2),hE_water
            tmp=tmp+1.0D0
        elseif(abs(Altitude_angle(2))>180.0D0) then
            write(*,'(A24)') "Return Flight"
            tmp=tmp+1.0D0
        elseif(CPU_time_min>20) then
            write(*,'(A24)') "Run time over"
            tmp=tmp+1.0D0
        elseif(isnan(hE)) then
            write(*,'(A24)') "Not A Numbar"
            tmp=tmp+1.0D0
        end if
        if(tmp>0.0D0) then
            write(*,'(A24,F12.5)') "Distance",sqrt(xE**2+yE**2)
            write(*,'(A24,I12)') "iteration",iteration
            write(*,'(A24,I12,A12,F12.3,A12)') "Computation Time",CPU_time_min,"[min]",CPU_time_s,"[s]"
            write(*,*) "------------------------------------------------------------"
            !write(*,*)
            iteration_max=iteration
            exit
        end if

        !偏差の総和を計算する
        sum_error=sum_error+(abs(PID_lon1(0))+abs(PID_lon2(0))+(abs(PID_lat1(0))+abs(PID_lat2(0)))*weighting_factor2)*time_step

        !入力の時間変化の総和を計算する
        sum_D_input=sum_D_input+(abs(Input(0)-Input_old(0))/de_max)
        sum_D_input=sum_D_input+(abs(Input(1)-Input_old(1))/dh_max)
        sum_D_input=sum_D_input+(abs(Input(2)-Input_old(2))/dr_max)

        !Input_oldを更新
        Input_old=Input

    end do

    !飛距離から偏差の総和を重みづけして差をとる
    !Distance=sqrt(Position(0)**2+Position(1)**2)-sum_error*weighting_factor1-sum_D_input*weighting_factor3
    Distance=sqrt(Position(0)**2+Position(1)**2)

end subroutine

subroutine Output_FlightLog(i_max_copy,iteration_max,FlightLog_data,FlightLog,genom_length,genom_list)
    !フライトログをcsvファイルで出力するための脳筋subroutine
    implicit none !暗黙の変数宣言を無効にする
    integer,intent(IN)::i_max_copy !実際に反復した回数
    integer,intent(IN)::iteration_max !最大反復回数
    integer,intent(IN)::FlightLog_data !FLightLogのデータ数
    real(8),dimension(0:iteration_max,0:FlightLog_data),intent(IN)::FlightLog !飛行解析のデータ
    integer,intent(IN)::genom_length !遺伝子情報の長さ
    real(8),dimension(0:genom_length-1),intent(IN)::genom_list !遺伝子情報
    integer::i

    open (10, file='FlightLog.csv', status='replace')
    do i = 0, i_max_copy
        write (10, '(F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,","&
        ,F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,","&
        ,F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,","&
        ,F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,","&
        ,F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3,",",F12.3)') &
        FlightLog(i,0),FlightLog(i,1),FlightLog(i,2),FlightLog(i,3),FlightLog(i,4)&
        ,FlightLog(i,5),FlightLog(i,6),FlightLog(i,7),FlightLog(i,8),FlightLog(i,9)&
        ,FlightLog(i,10),FlightLog(i,11),FlightLog(i,12),FlightLog(i,13),FlightLog(i,14)&
        ,FlightLog(i,15),FlightLog(i,16),FlightLog(i,17),FlightLog(i,18),FlightLog(i,19)&
        ,FlightLog(i,20),FlightLog(i,21),FlightLog(i,22),FlightLog(i,23),FlightLog(i,24)&
        ,FlightLog(i,25),FlightLog(i,26),FlightLog(i,27),FlightLog(i,28),FlightLog(i,29)&
        ,FlightLog(i,30),FlightLog(i,31),FlightLog(i,32),FlightLog(i,33),FlightLog(i,34)&
        ,FlightLog(i,35),FlightLog(i,36),FlightLog(i,37),FlightLog(i,38),FlightLog(i,39)&
        ,FlightLog(i,40),FlightLog(i,41),FlightLog(i,42),FlightLog(i,43),FlightLog(i,44)&
        ,FlightLog(i,45),FlightLog(i,46),FlightLog(i,47)
    end do
    close (10)

    open (10, file='List.csv', status='replace')
    write(10,*) genom_list(0),",",genom_list(1),",",genom_list(2)
    do i=3,8
        write(10,*) genom_list(i) !1行ずつ値を読み込む
    end do
    do i=0,3
        write(10,*) genom_list(9+i*3),",",genom_list(10+i*3),",",genom_list(11+i*3)
    end do
    close(10) !ファイルを閉じる

end subroutine

subroutine GA_NLFS
    !遺伝的アルゴリズムによって飛行経路を最適化する
    !符号:実数値エンコーディング
    !選択:エリート主義
    !交叉:二点交叉(or 一様交叉)
    !変位:摂動
    !世代:定常状態モデル
    implicit none !暗黙の変数宣言を無効にする
    integer,parameter::genom_length = 21 !遺伝子情報の長さ
    integer,parameter::max_genom_list = 50 !遺伝子集団の大きさ
    integer,parameter::select_genom = 20 !遺伝子選択数
    integer,parameter::max_generation = 50 !繰り返す世代
    real(8),parameter::individual_mutation = 0.010D0 !固体突然変異率
    real(8),parameter::genom_mutation = 0.010D0 !遺伝子突然変異率
    integer,parameter::FlightLog_data = 47 !飛行解析のデータ数
    real(8),dimension(0:genom_length-1)::genom_list !遺伝子情報
    real(8),dimension(0:genom_length-1,0:max_genom_list-1)::genom_class !遺伝子集団
    real(8),dimension(0:genom_length-1,0:select_genom-1)::elite_class !エリート
    real(8),dimension(0:genom_length-1,0:select_genom-1)::progeny_class !子孫
    real(8),dimension(0:max_genom_list-1)::Distance !飛行距離 [m]
    real(8)::time_max !最大飛行時間 [s]
    real(8)::time_step !タイムステップ [s]
    integer::iteration_max !最大反復回数 [s]
    integer::i_max_copy !最大反復回数のコピー [s]
    integer::pilot_method !操縦方式．0：エレベータ，1：重心移動
    !カウンター
    integer::i,j,generation
    !計算時間計測
    real(8)::CPU_time_s
    integer::initial_time,finish_time,time_rate,CPU_time_min,CPU_time_h !initial time, finish time, time rate
    
    call system_clock(initial_time) !開始時間の読み込み

    write(*,*) "------------------------------------------------------------"
    write(*,*)
    write(*,'(A24,I12,A12,I12)') "Generaion",0,"/",max_generation 
    call system_clock(finish_time,time_rate) !終了時間の読み込み
    CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
    if(CPU_time_s>60) then
        CPU_time_min=int(CPU_time_s/60.0D0)
        CPU_time_s=CPU_time_s-60.D0*CPU_time_min
    end if
    write(*,'(A24,I12,A12,F12.3,A12)') "Computation Time",CPU_time_min,"[min]",CPU_time_s,"[s]"
    write(*,*) "------------------------------------------------------------"

    !------------------------------------------------------!
    !--------------------遺伝子集団を作る--------------------!
    !------------------------------------------------------!

    do i=0,max_genom_list-1
        call create_genom(genom_length,genom_list)
        do j=0,genom_length-1
            genom_class(j,i)=genom_list(j)
        end do
    end do
    
    !------------------------------------------------------!
    !--------------------遺伝子を評価する--------------------!
    !------------------------------------------------------!

    call evalution(genom_length,max_genom_list,genom_class,Distance,FlightLog_data)

    !-------------------------------------------------------!
    !--------------------遺伝子の並び替え--------------------!
    !-------------------------------------------------------!

    call sort_genom(max_genom_list,genom_length,Distance,genom_class)

    !---------------------------------------------------!
    !--------------------世代のループ--------------------!
    !---------------------------------------------------!

    do generation=0,max_generation

        write(*,*) "------------------------------------------------------------"
        write(*,'(A24,I12,A12,I12)') "Generaion",generation,"/",max_generation 
        call system_clock(finish_time,time_rate) !終了時間の読み込み
        CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
        if(CPU_time_s>60) then
            CPU_time_min=int(CPU_time_s/60.0D0)
            CPU_time_s=CPU_time_s-60.D0*CPU_time_min
        end if
        write(*,'(A24,I12,A12,F12.3,A12)') "Computation Time",CPU_time_min,"[min]",CPU_time_s,"[s]"
        write(*,*) "------------------------------------------------------------"

        !-----------------------------------------------------!
        !--------------------エリートの選択--------------------!
        !-----------------------------------------------------!

        call select_elite(select_genom,max_genom_list,genom_length,genom_class,elite_class)

        !-------------------------------------------!
        !--------------------交叉--------------------!
        !-------------------------------------------!

        call cross_over(select_genom,genom_length,max_genom_list,genom_class,elite_class,progeny_class)
        
        !-----------------------------------------------!
        !--------------------突然変異--------------------!
        !-----------------------------------------------!

        call mutation(individual_mutation,genom_mutation,genom_length,max_genom_list,genom_class)

        !------------------------------------------------------!
        !--------------------遺伝子を評価する--------------------!
        !------------------------------------------------------!

        call evalution(genom_length,max_genom_list,genom_class,Distance,FlightLog_data)

        !-------------------------------------------------------!
        !--------------------遺伝子の並び替え--------------------!
        !-------------------------------------------------------!

        call sort_genom(max_genom_list,genom_length,Distance,genom_class)

        !-------------------------------------------------!
        !--------------------結果を表示--------------------!
        !-------------------------------------------------!

        if(mod(generation,1)==0 .OR. generation==max_generation)  then
            call output_genom(genom_length,max_genom_list,genom_class,Distance,generation,FlightLog_data)
        end if

    end do

    write(*,*) "------------------------------------------------------------"
    call system_clock(finish_time,time_rate) !終了時間の読み込み
    CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
    if(CPU_time_s>60) then
        CPU_time_min=int(CPU_time_s/60.0D0)
        CPU_time_s=CPU_time_s-60.D0*CPU_time_min
        if(CPU_time_min>60) then
            CPU_time_h=int(CPU_time_min/60)
            CPU_time_min=CPU_time_min-60*CPU_time_h
        end if
    end if
    write(*,'(A24,I12,A12,I12,A12,F12.3,A12)') "Computation Time",CPU_time_h,"[h]",CPU_time_min,"[min]",CPU_time_s,"[s]"

end subroutine

subroutine create_genom(genom_length,genom_list)
    !ランダムに遺伝子を生成する
    implicit none !暗黙の変数宣言を無効にする
    !引数
    integer,intent(IN)::genom_length !遺伝子情報の長さ
    real(8),dimension(0:genom_length-1),intent(OUT)::genom_list !遺伝子情報
    !変数
    real(8),parameter::phi_max =0.0D0 !最大初期バンク角 [deg]
    real(8),parameter::theta_max =0.0D0 !最大初期ピッチ角 [deg]
    real(8),parameter::psi_max =0.0D0 !最大初期ヨー角 [deg]
    real(8)::time_max !最大飛行時間 [s]
    real(8),parameter::max_dive_angle =-20.0D0 !最大ダイブ角 [deg]
    real(8),parameter::glide_angle =-2.000D0 !定常滑空角 [deg]
    real(8),parameter::max_bank_angle =0.000D0 !最大バンク角 [deg]
    real(8)::rnd !0~1の乱数
    real(8)::tmp
    !カウンター
    integer::i,j,k

    !操縦方式，飛行時間，タイムステップの読み込み
    open(10,file="LFS.txt") !ファイルを開く
    read(10,*)  !操縦方式
    read(10,*) time_max !最大飛行時間
    read(10,*)  !タイムステップ
    close(10) !ファイルを閉じる
    
    call random(rnd)
    genom_list(0)=phi_max*rnd !初期バンク角 [deg]
    call random(rnd)
    if(rnd<0.5D0) then !半分の確率で負のバンク角にする
        genom_list(0)=-genom_list(0)
    end if

    call random(rnd)
    genom_list(1)=theta_max*rnd !初期ピッチ角 [deg]
    call random(rnd)
    if(rnd<0.5D0) then !半分の確率で負のピッチ角にする
        genom_list(1)=-genom_list(1)
    end if

    call random(rnd)
    genom_list(2)=psi_max*rnd !初期ヨー角 [deg]
    call random(rnd)
    if(rnd<0.5D0) then !半分の確率で負のヨー角にする
        genom_list(2)=-genom_list(2)
    end if

    call random(rnd)
    genom_list(3)=time_max*rnd*0.25D0 !引き起こし時刻 [s]

    call random(rnd)
    genom_list(4)=time_max*rnd*0.5D0 !バンク開始時刻 [s]

    call random(rnd)
    genom_list(5)=time_max*rnd !バンク終了時刻 [s]
    if(genom_list(5)<genom_list(4)) then !開始時刻＜終了時刻となるように並び替え
        tmp=genom_list(4)
        genom_list(4)=genom_list(5)
        genom_list(5)=tmp
    end if

    call random(rnd)
    genom_list(6)=glide_angle+max_dive_angle*rnd !ダイブ角 [deg]
    !genom_list(6)=max_dive_angle*rnd !ダイブ角 [deg]

    call random(rnd)
    !genom_list(7)=glide_angle !定常滑空角 [deg]
    !genom_list(7)=max_dive_angle*rnd !定常滑空角 [deg]
    genom_list(7)=glide_angle*rnd !定常滑空角 [deg]

    call random(rnd)
    genom_list(8)=max_bank_angle*rnd !バンク角 [deg]
    call random(rnd)
    if(rnd<0.5D0) then !半分の確率で負のバンク角にする
        genom_list(8)=-genom_list(8)
    end if

    do j=0,11
        call random(rnd)
        genom_list(9+j)=rnd !PIDパラメータ
    end do

end subroutine

subroutine evalution(genom_length,max_genom_list,genom_class,Distance,FlightLog_data)
    !各遺伝子を評価する．
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan !プログラム中で使用する関数の名前を宣言
    !引数
    integer,intent(IN)::genom_length !遺伝子情報の長さ
    integer,intent(IN)::max_genom_list !遺伝子集団の大きさ
    real(8),dimension(0:genom_length-1,0:max_genom_list-1),intent(IN)::genom_class !遺伝子集団
    real(8),dimension(0:max_genom_list-1),intent(OUT)::Distance !飛行距離 [m]
    integer,intent(IN)::FlightLog_data !FLightLogのデータ数
    !変数
    real(8),dimension(0:genom_length-1)::genom_list !遺伝子情報
    real(8),dimension(:,:),allocatable::FlightLog !飛行解析のデータ
    real(8)::time_max !最大飛行時間 [s]
    real(8)::time_step !タイムステップ [s]
    integer::iteration_max !最大反復回数 [s]
    integer::i_max_copy !最大反復回数のコピー [s]
    integer::pilot_method !操縦方式．0：エレベータ，1：重心移動
    !カウンター
    integer::i,j

    !飛行時間，タイムステップの読み込み
    open(10,file="LFS.txt") !ファイルを開く
    read(10,*) pilot_method !操縦方式
    read(10,*) time_max !最大飛行時間
    read(10,*) time_step !タイムステップ
    close(10) !ファイルを閉じる

    iteration_max=int(time_max/time_step)+1 !最大反復回数の読み込み
    allocate(FlightLog(0:iteration_max,0:FlightLog_data))

    do i=0,max_genom_list-1

        !遺伝子情報の読み込み
        do j=0,genom_length-1
            genom_list(j)=genom_class(j,i) !1行ずつ値を読み込む
        end do

        iteration_max=int(time_max/time_step)+1 !最大反復回数の読み込み
        i_max_copy=iteration_max
        write(*,'(A24,I12,A12,I12)') "Genom",i+1,"/",max_genom_list
        call NLFS(Distance(i),genom_length,genom_list,i_max_copy,pilot_method,FlightLog_data,FlightLog)
        if(isnan(Distance(i))) stop
        !write(*,'(A12,F12.5)') "Distance",Distance(i)
    
    end do

end subroutine

subroutine sort_genom(max_genom_list,genom_length,Distance,genom_class)
    !ゲノムを優秀な順番に並び替える
    implicit none
    integer,intent(IN):: max_genom_list
    integer,intent(IN):: genom_length
    real(8),dimension(0:genom_length-1,0:max_genom_list-1),intent(INOUT)::genom_class !遺伝子集団
    real(8),dimension(0:max_genom_list-1),intent(INOUT):: Distance
    real(8),dimension(0:genom_length-1)::genom_list !遺伝子情報
    real(8)::tmp
    integer::i,j,k

    !0番目の要素から
    do i=1,max_genom_list-1
        do j=i,1,-1
            if(Distance(j)>Distance(j-1)) then
                !Distanceを並び替え
                tmp=Distance(j-1)
                Distance(j-1)=Distance(j)
                Distance(j)=tmp
                !genom_classを並び替え
                do k=0,genom_length-1
                    genom_list(k)=genom_class(k,j-1)
                    genom_class(k,j-1)=genom_class(k,j)
                    genom_class(k,j)=genom_list(k)
                end do
            end if
        end do
    end do

end subroutine

subroutine select_elite(sselect_genom,max_genom_list,genom_length,genom_class,elite_class)
    !エリートを選択する
    implicit none
    integer,intent(IN)::sselect_genom !遺伝子選択数
    integer,intent(IN):: max_genom_list
    integer,intent(IN):: genom_length
    real(8),dimension(0:genom_length-1,0:max_genom_list-1),intent(INOUT)::genom_class !遺伝子集団
    real(8),dimension(0:genom_length-1,0:sselect_genom-1),intent(OUT)::elite_class !エリート
    real(8)::tmp
    integer::i,j,k

    do i=0,sselect_genom-1
        do j=0,genom_length-1
            elite_class(j,i)=genom_class(j,i)
        end do
    end do

end subroutine

subroutine cross_over(sselect_genom,genom_length,max_genom_list,genom_class,elite_class,progeny_class)
    !ランダムに2つのエリートを選択し，一様交叉を行う
    implicit none
    integer,intent(IN)::sselect_genom !遺伝子選択数
    integer,intent(IN)::genom_length
    integer,intent(IN)::max_genom_list
    real(8),dimension(0:genom_length-1,0:max_genom_list-1),intent(INOUT)::genom_class !遺伝子集団
    real(8),dimension(0:genom_length-1,0:sselect_genom-1),intent(IN)::elite_class !エリート
    real(8),dimension(0:genom_length-1,0:sselect_genom-1),intent(OUT)::progeny_class !子孫
    real(8),dimension(0:genom_length-1)::genom_list1,genom_list2 !遺伝子情報
    integer::tmp,genom1,genom2
    real(8)::rnd
    integer::i,j,k

    !select_genomの数だけ子孫を作る
    do i=0,sselect_genom-1
        !ランダムにエリートを2つ選択する
        call random(rnd)
        do j=0,genom_length-1
            genom_list1(j)=elite_class(j,int(sselect_genom*rnd))
        end do
        call random(rnd)
        do j=0,genom_length-1
            genom_list2(j)=elite_class(j,int(sselect_genom*rnd))
        end do
        !一様交叉を行う
        !do j=0,genom_length-1
        !    call random(rnd)
        !    If(rnd<0.5D0) then
        !        progeny_class(j,i)=genom_list1(j)
        !    else
        !        progeny_class(j,i)=genom_list2(j)
        !    end if
        !二点交叉を行う
        !ランダムに2つの遺伝子を選択する
        call random(rnd)
        genom1=int(genom_length*rnd)
        call random(rnd)
        genom2=int(genom_length*rnd)
        if(genom1>genom2) then !genom1<genom2になるよう並び替え
            tmp=genom1
            genom1=genom2
            genom2=tmp
        end if
        do j=0,genom_length-1 !まずgenom_list1を格納し
            progeny_class(j,i)=genom_list1(j)
        end do
        do j=genom1,genom2 !指定した範囲だけgenom_list2と入れ替える
            progeny_class(j,i)=genom_list2(j)
        end do 
    end do

    !現在の世代の評価の低い個体と子孫を入れ替える
    do i=0,sselect_genom-1
        do j=0,genom_length-1
            genom_class(j,max_genom_list-1-i)=progeny_class(j,i)
        end do
    end do

end subroutine

subroutine mutation(individual_mutation,genom_mutation,genom_length,max_genom_list,genom_class)
    !突然変異を行う
    implicit none
    real(8),intent(IN)::individual_mutation !固体突然変異率
    real(8),intent(IN)::genom_mutation !遺伝子突然変異率
    integer,intent(IN)::genom_length
    integer,intent(IN)::max_genom_list
    real(8),dimension(0:genom_length-1,0:max_genom_list-1),intent(INOUT)::genom_class !遺伝子集団
    real(8),dimension(0:genom_length-1)::genom_list !遺伝子情報
    real(8)::time_max !最大飛行時間 [s]
    real(8),parameter::max_dive_angle =-60.0D0 !最大ダイブ角 [deg]
    real(8),parameter::max_glide_angle =-2.0D0 !最大定常滑空角 [deg]
    real(8),parameter::max_bank_angle =10.0D0 !最大バンク角 [deg]
    real(8)::rnd
    real(8)::tmp
    integer::i,j,k

    !操縦方式，飛行時間，タイムステップの読み込み
    open(10,file="LFS.txt") !ファイルを開く
    read(10,*)  !操縦方式
    read(10,*) time_max !最大飛行時間
    read(10,*)  !タイムステップ
    close(10) !ファイルを閉じる

    !個体のループ
    do i=0,max_genom_list-1
        call random(rnd)
        if(rnd<individual_mutation) then !個体突然変異
            call create_genom(genom_length,genom_list)
            do j=0,genom_length-1
                genom_class(j,i)=genom_list(j)
            end do
        end if
        
        !遺伝子のループ
        call create_genom(genom_length,genom_list)
        do j=0,genom_length-1
            call random(rnd)
            if(rnd<genom_mutation) then !遺伝子突然変異率
                genom_class(j,i)=genom_list(j)
            end if
        end do
        
        !バンク開始時刻＜終了時刻となるように並び替え
        if(genom_class(5,i)<genom_class(4,i)) then 
            tmp=genom_Class(4,i)
            genom_class(4,i)=genom_class(4,i)
            genom_class(5,i)=tmp
        end if
            
    end do

end subroutine

subroutine output_genom(genom_length,max_genom_list,genom_class,Distance,generation,FlightLog_data)
    !結果の出力
    implicit none
    integer,intent(IN)::genom_length
    integer,intent(IN)::max_genom_list
    real(8),dimension(0:genom_length-1,0:max_genom_list-1),intent(IN)::genom_class !遺伝子集団
    real(8),dimension(0:max_genom_list-1),intent(IN)::Distance !飛行距離 [m]
    integer,intent(IN)::generation
    integer,intent(IN)::FlightLog_data !FLightLogのデータ数
    real(8),dimension(0:genom_length-1)::genom_list !遺伝子情報
    real(8)::Distance_Max,Distance_Min,Distance_Ave
    real(8),dimension(:,:),allocatable::FlightLog !飛行解析のデータ
    real(8)::time_max !最大飛行時間 [s]
    real(8)::time_step !タイムステップ [s]
    integer::iteration_max !最大反復回数 [s]
    integer::i_max_copy !最大反復回数のコピー [s]
    integer::pilot_method !操縦方式．0：エレベータ，1：重心移動
    real(8)::rnd,tmp
    integer::i,j,k

    Distance_Max=Distance(0)
    Distance_Min=Distance(max_genom_list-1)
    do i=0,max_genom_list-1
        Distance_Ave=Distance_Ave+Distance(i)
    end do
    Distance_Ave=Distance_Ave/max_genom_list
    do i=0,genom_length-1
        genom_list(i)=genom_class(i,0)
    end do

    !飛行時間，タイムステップの読み込み
    open(10,file="NLFS.txt") !ファイルを開く
    read(10,*) pilot_method !操縦方式
    read(10,*) time_max !最大飛行時間
    read(10,*) time_step !タイムステップ
    close(10) !ファイルを閉じる

    iteration_max=int(time_max/time_step)+1 !最大反復回数の読み込み
    allocate(FlightLog(0:iteration_max,0:FlightLog_data))

    i_max_copy=iteration_max
    call NLFS(tmp,genom_length,genom_list,iteration_max,pilot_method,FlightLog_data,FlightLog)

    call Output_FlightLog(i_max_copy,iteration_max,FlightLog_data,FlightLog,genom_length,genom_list) !飛行解析データをcsvファイルで書き出し

    write(*,*) "-----------------------------------------------------------------------------------"
    write(*,'(A15,A1,I3)') "generation","=",generation
    write(*,'(A15,A1,F8.3)') "MAX","=",Distance_Max
    write(*,'(A15,A1,F8.3)') "MIN","=",Distance_Min
    write(*,'(A15,A1,F8.3)') "AVERAGE","=",Distance_Ave
    write(*,'(A15)') "BEST GENOM"
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "phi","=",genom_list(0),"theta","=",genom_list(1),"psi","=",genom_list(2)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "time_pullup","=",genom_list(3),"dive_angle","=",genom_list(6),"cluse_angle","=",genom_list(7)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "time_bank","=",genom_list(4),"time_bank_end","=",genom_list(5),"bank_angle","=",genom_list(8)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "PID_lon1[0]","=",genom_list(9),"PID_lon2[0]","=",genom_list(12)&
    ,"PID_lat1[0]","=",genom_list(15),"PID_lat2[0]","=",genom_list(18)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "PID_lon1[1]","=",genom_list(10),"PID_lon2[1]","=",genom_list(13)&
    ,"PID_lat1[1]","=",genom_list(16),"PID_lat2[1]","=",genom_list(19)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "PID_lon1[1]","=",genom_list(11),"PID_lon2[2]","=",genom_list(14)&
    ,"PID_lat1[2]","=",genom_list(17),"PID_lat2[2]","=",genom_list(20)

end subroutine

subroutine random(rnd)
    !疑似的な乱数を生成する
    !https://qiita.com/ocian/items/e5eabe60c6b31cb48fac
    implicit none
    integer :: seedsize
    integer,allocatable :: seed(:)
    real(8),intent(OUT) :: rnd
    integer :: c !時間を入れる
    integer :: i !カウンター

    call system_clock(count=c) !時間を取得
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    call random_seed(get=seed)
    seed = c !時間を全部に代入
    call random_number(rnd) !rndに乱数をセット
end subroutine

subroutine Plane_calculation &
        (Altitude_angle,Ground_speed,Angular_velocity,gust,Input &
        ,hE,Lift,Drag,Force,Moment,epsilon,Acceleration,Angular_acceleration)
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan,sqrt,sin,cos,tan !プログラム中で使用する関数の名前を宣言
    !パラメータ
    real(8),parameter::g =9.795 !重力加速度 [m/s^2]
    !引数
    real(8),dimension(0:2),intent(IN)::Altitude_angle !姿勢角 [Φ,θ,Ψ] [deg]
    real(8),dimension(0:2),intent(INOUT)::Ground_speed !機体軸における対地速度ベクトル [u,v,w] [m/s]
    real(8),dimension(0:2),intent(INOUT)::Angular_velocity !機体軸における角速度ベクトル [p,q,r] [deg/s]
    real(8),dimension(0:2),intent(IN)::Input !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg] [δe,δh,δr]
    real(8),dimension(0:2),intent(IN)::gust !機体軸における突風ベクトル [ug,vg,wg] [m/s]
    real(8),dimension(0:2),intent(OUT)::Acceleration !機体軸における対地加速度ベクトル [du/dt,dv/dt,dw/dt] [m/s]
    real(8),dimension(0:2),intent(OUT)::Angular_acceleration !機体軸における角加速度ベクトル [dp/dt,dq/dt,dr/dt] [deg/s]
    real(8)::hE !高度 [m]
    real(8),intent(OUT)::Lift !揚力 [N]
    real(8),intent(OUT)::Drag !抗力 [N]
    real(8),dimension(0:2),intent(OUT)::Force !機体にはたらく外力 [X,Y,Z] [N]
    real(8),dimension(0:2),intent(OUT)::Moment !機体にはたらくモーメント [L,M,N] [N*m]
    real(8),intent(INOUT)::epsilon !主翼位置での吹きおろし角 [deg]
    !機体データ
    integer::span_div,chord_div !スパン方向分割数
    real(8)::Mass !機体重量 [kg]
    real(8)::Ixx,Iyy,Izz,Ixz !慣性モーメント [kg*m^4]
    real(8)::CD_fuselage,S_fuselage,chord_fuselage !カウル，胴体の抗力係数に(1/2)*ρ*Sを掛け合わせたもの [kg/m]
    real(8),dimension(0:14)::Cdp_coef !揚力係数，抗力係数，ピッチングモーメント係数を多項式近似したときの係数
    !その他
    real(8)::Drag_fuselage,CD0
    real(8)::rho,mu,Re
    real(8)::Velocity
    real(8)::X,Y,Z,L,M,N
    real(8)::phi,theta,psi
    real(8)::de,dh,dr
    real(8)::p,q,r
    real(8)::u,v,w,ug,vg,wg
    real(8)::alpha,beta,gamma
    real(8)::pi,deg,rad !円周率，degへの変換，radへの変換
    real(8)::tmp
    !カウンター
    integer::i

    !write(*,'(A24)') "Plane_calculation"

    pi=4.0D0*atan(1.0D0) !円周率
    deg=(180.0D0/pi) !degへの変換
    rad=(pi/180.0D0) !radへの変換

    Lift=0.0D0 !揚力 [N]
    Drag=0.0D0 !抗力 [N]
    Force=0.0D0 !機体にはたらく外力 [X,Y,Z] [N]
    Moment=0.0D0 !機体にはたらくモーメント [L,M,N] [N*m]
    Acceleration=0.0D0
    Angular_acceleration=0.0D0

    !機体データ読み込み
    open(10,file="Plane_data.txt") !ファイルを開く
        read(10,*) Mass !機体重量 [kg]
        read(10,*) Ixx,Iyy,Izz,Ixz !慣性モーメント [kg*m^4]
        read(10,*) CD_fuselage,S_fuselage,chord_fuselage !カウル，胴体の抗力係数に(1/2)*ρ*Sを掛け合わせたもの [kg/m]
    close(10) !ファイルを閉じる

    !わかりやすいように変数変換
    u=Ground_speed(0);          w=Ground_speed(2);              v=Ground_speed(1)           !u,v,w [m/s]
    ug=gust(0);                 wg=gust(2);                     vg=gust(1)                  !u,v,w [m/s]
    phi=Altitude_angle(0);      theta=Altitude_angle(1);            psi=Altitude_angle(2)   !Φ，θ，ψ [rad]を格納
    p=Angular_velocity(0);      q=Angular_velocity(1);          r=Angular_velocity(2)       !p,q,r [deg]       
    de=Input(0);                dh=Input(1);                    dr=Input(2)                 !操舵を格納

    !速度V [m/s]，全機迎角α [deg]，横滑り角β [deg]を計算する
    Velocity=sqrt((u+ug)**2+(v+vg)**2+(w+wg)**2)
    alpha=deg*asin((w+wg)/(u+ug))
    beta=deg*asin((v+vg)/Velocity)
    gamma=theta-alpha

    !水平尾翼計算
    open(10,file="Tail_data.txt") !ファイルを開く
        read(10,*) span_div,chord_div !スパン方向分割数読み込み
    close(10) !ファイルを閉じる
    call Tail_calculation(Lift,Drag,Force,Moment,Velocity,alpha,q,de,dh,epsilon,span_div,chord_div)
    !if(isnan(Lift)) write(*,'(A24)') "Tail"

    !if(1==0) then
        !write(*,'(A24)') "Tail_calculation"
        !write(*,'(A12,F12.3)') "hE",hE
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "u",u,"v",v,"w",w
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "phi",phi,"theta",theta,"psi",psi
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "p",p,"q",q,"r",r
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "alpha",alpha,"beta",beta,"gamma",gamma
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "de",de,"dh",dh,"dr",dr
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "X",Force(0),"Y",Force(1),"Z",Force(2)
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "L",Moment(0),"M",Moment(1),"N",Moment(2)
        !write(*,'(A12,F12.5,A12,F12.5,A12,F12.5)')  &
        !"du/dt",Acceleration(0),"dv/dt",Acceleration(1),"dw/dt",Acceleration(2)
        !write(*,'(A12,F12.5,A12,F12.5,A12,F12.5)')  &
        !"dp/dt",Angular_Acceleration(0),"dq/dt",Angular_Acceleration(1),"dr/dt",Angular_Acceleration(2)
        !write(*,*)
    !end if

    !主翼計算    
    open(10,file="Wing_data.txt") !ファイルを開く
        read(10,*) span_div,tmp !スパン方向分割数読み込み
    close(10) !ファイルを閉じる
    !call Wing_calculation(Lift,Drag,Force,Moment,Velocity,alpha,beta,p,r,hE,dh,epsilon,span_div)
    call Wing_calculation_linear(Lift,Drag,Force,Moment,Velocity,alpha,beta,p,r,hE,dh,epsilon,span_div)
    !if(isnan(Lift)) write(*,'(A24)') "Wing"

    !if(1==0) then
        !write(*,'(A24)') "Wing_calculation"
        !write(*,'(A12,F12.3)') "hE",hE
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "u",u,"v",v,"w",w
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "phi",phi,"theta",theta,"psi",psi
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "p",p,"q",q,"r",r
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "alpha",alpha,"beta",beta,"gamma",gamma
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "de",de,"dh",dh,"dr",dr
        !write(*,'(A12,F12.3,A12,F12.3)') "Lift",Lift,"Drag",Drag
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "X",Force(0),"Y",Force(1),"Z",Force(2)
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "L",Moment(0),"M",Moment(1),"N",Moment(2)
        !write(*,'(A12,F12.5,A12,F12.5,A12,F12.5)')  &
        !"du/dt",Acceleration(0),"dv/dt",Acceleration(1),"dw/dt",Acceleration(2)
        !write(*,'(A12,F12.5,A12,F12.5,A12,F12.5)')  &
        !"dp/dt",Angular_Acceleration(0),"dq/dt",Angular_Acceleration(1),"dr/dt",Angular_Acceleration(2)
        !write(*,*)
    !end if

    !垂直尾翼計算
    open(10,file="Fin_data.txt") !ファイルを開く
        read(10,*) span_div,chord_div !スパン方向分割数読み込み
        read(10,*) rho,mu !空気密度 [kg/m^3]，動粘性係数 [m^2/s]
        do i=0,5
            read(10,*) tmp
        end do
        read(10,*) Cdp_coef(0:14) !抗力係数を多項式近似したときの係数
    close(10) !ファイルを閉じる
    call Fin_calculation(Lift,Drag,Force,Moment,Velocity,beta,p,r,dr,dh,span_div,chord_div)

    !if(1==0) then
        !write(*,'(A24)') "Fin_calculation"
        !write(*,'(A12,F12.3)') "hE",hE
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "u",u,"v",v,"w",w
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "phi",phi,"theta",theta,"psi",psi
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "p",p,"q",q,"r",r
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "alpha",alpha,"beta",beta,"gamma",gamma
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "de",de,"dh",dh,"dr",dr
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "X",Force(0),"Y",Force(1),"Z",Force(2)
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "L",Moment(0),"M",Moment(1),"N",Moment(2)
        !write(*,'(A12,F12.5,A12,F12.5,A12,F12.5)')  &
        !"du/dt",Acceleration(0),"dv/dt",Acceleration(1),"dw/dt",Acceleration(2)
        !write(*,'(A12,F12.5,A12,F12.5,A12,F12.5)')  &
        !"dp/dt",Angular_Acceleration(0),"dq/dt",Angular_Acceleration(1),"dr/dt",Angular_Acceleration(2)
        !write(*,*)
    !end if

    !カウル抗力計算
    Re=(Velocity*1.0D0)/mu
    !垂直尾翼の迎角0 [deg]での抗力係数を計算
    CD0=Cdp_coef(0)
    do i=9,14
        tmp=((Re/100000.0D0)**(i-8))
        CD0=CD0+Cdp_coef(i)*tmp
    end do
    CD0=CD_fuselage/CD0 !垂直尾翼の抗力係数の係数をカウルに適用できるようにする
    !カウルの抗力係数を計算
    CD_fuselage=0.0D0
    do i=0,8
        tmp=alpha**i
        CD_fuselage=CD_fuselage+CD0*Cdp_coef(i)*tmp
    end do
    do i=9,14
        tmp=((Re/100000.0D0)**(i-8))
        CD_fuselage=CD_fuselage+CD0*Cdp_coef(i)*tmp
    end do
    Drag_fuselage=0.5D0*rho*Velocity**2*CD_fuselage*S_fuselage
    Drag=Drag+Drag_fuselage
    !Force(0)=Force(0)-Drag_fuselage*cos(rad*alpha)
    !Force(2)=Force(2)+Drag_fuselage*sin(rad*alpha)
    Force(0)=Lift*sin(rad*alpha)-Drag*cos(rad*alpha)
    Force(2)=-Lift*cos(rad*alpha)+Drag*sin(rad*alpha)

    X=Force(0)
    Y=Force(1)
    Z=Force(2)
    L=(Iyy-Izz)*rad*q*r+Ixz*rad*p*q+deg*Moment(0)
    M=(Izz-Ixx)*rad*r*p+Ixz*rad*(r*r-p*p)+deg*Moment(1)
    N=(Ixx-Iyy)*rad*p*q-Ixz*rad*q*r+deg*Moment(2)

    !非線形6自由度運動方程式を解く．青本p.16
    Acceleration(0)=-rad*q*w+rad*r*v-g*sin(rad*theta)+(X/Mass)
    Acceleration(1)=-rad*r*u+rad*p*w+g*cos(rad*theta)*sin(rad*phi)+(Y/Mass)
    Acceleration(2)=-rad*p*v+rad*q*u+g*cos(rad*theta)*cos(rad*phi)+(Z/Mass)
    !write(*,"(A24,F12.5)") "Mass",Mass
    !write(*,"(A24,F12.5,A24,F12.5)") "(Z/Mass)",(Z/Mass),"g*cos(rad*theta)*cos(rad*phi)",g*cos(rad*theta)*cos(rad*phi)
    Angular_acceleration(0)=((L/Ixx)+(Ixz/Ixx)*(N/Izz))/(1.0D0-(Ixz**2)/(Izz*Ixx))
    Angular_acceleration(1)=M/Iyy
    Angular_acceleration(2)=((N/Izz)+(Ixz/Izz)*(L/Ixx))/(1.0D0-(Ixz**2)/(Ixx*Izz))

    !if(1==0) then
        !write(*,*)
        !write(*,'(A12,F12.3)') "hE",hE
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "u",u,"v",v,"w",w
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "phi",phi,"theta",theta,"psi",psi
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "p",p,"q",q,"r",r
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "alpha",alpha,"beta",beta,"gamma",gamma
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "de",de,"dh",dh,"dr",dr
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "X",Force(0),"Y",Force(1),"Z",Force(2)
        !write(*,'(A12,F12.3,A12,F12.3,A12,F12.3)') "L",Moment(0),"M",Moment(1),"N",Moment(2)
        !write(*,'(A12,F12.5,A12,F12.5,A12,F12.5)')  &
        !"du/dt",Acceleration(0),"dv/dt",Acceleration(1),"dw/dt",Acceleration(2)
        !write(*,'(A12,F12.5,A12,F12.5,A12,F12.5)')  &
        !"dp/dt",Angular_Acceleration(0),"dq/dt",Angular_Acceleration(1),"dr/dt",Angular_Acceleration(2)
        !write(*,*)
    !end if

end subroutine

subroutine Wing_calculation(Lift,Drag,Force,Moment,Velocity,alpha,beta,p,r,hE,dh,epsilon,span_div)
    !機体の姿勢角による翼素の高度への影響は無視．横滑りによって渦糸が斜めに出る影響は無視．
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan,sqrt,sin,cos,tan !プログラム中で使用する関数の名前を宣言
    !引数
    real(8),intent(INOUT)::Lift !揚力 [N]
    real(8),intent(INOUT)::Drag !抗力 [N]
    real(8),dimension(0:2),intent(INOUT)::Force !機体にはたらく外力 [X,Y,Z] [N]
    real(8),dimension(0:2),intent(INOUT)::Moment !機体にはたらくモーメント [L,M,N] [N*m]
    real(8),intent(IN)::Velocity,alpha,beta,p,r,hE,dh !対気速度 [m/s]，迎角 [deg]，横滑り角 [deg]，角速度 [deg/s]高度 [m]，重心移動距離 [-]
    real(8),intent(OUT)::epsilon !主翼位置での吹きおろし角 [deg]
    integer,intent(INOUT)::span_div
    !LLT
    !パラメータ
    real(8),parameter::circulation_coef =0.10D0 !循環の更新に使う係数
    real(8),parameter::error = 0.0005D0 !収束判定に使う誤差
    integer::iteration_max=2000 !最大反復計算回数
    real(8),parameter::alpha_max=20.0D0 !計算を発散させないためのalphaの最大値
    real(8),parameter::alpha_min=-10.0D0 !計算を発散させないためのalphaの最小値
    real(8),parameter::Re_max =1000000.0D0 !計算を発散させないためのRe数の最小値
    real(8),parameter::Re_min =100000.0D0 !計算を発散させないためのRe数の最小値
    !環境変数
    real(8)::rho,mu !空気密度 [kg/m^3]，動粘性係数 [m^2/s]
    !解析用
    real(8)::chord_div !スパン方向分割数，コード方向分割数
    !翼の幾何学形状
    real(8)::dy,ds !翼素間隔，パネルの反復 [m]
    real(8),dimension(0:2*span_div-1)::dz !高さの差 [m]
    real(8)::hac,hspar !空力中心位置，桁位置
    real(8)::hcg !重心位置
    real(8)::chord_mac !平均空力翼弦長 [m]
    integer::num_mac !平均空力翼弦長になるリブ番号
    real(8),dimension(0:8)::zc !中心線を多項式近似したときの係数
    real(8),dimension(0:2*span_div)::chord !コード長 [m]
    real(8),dimension(0:2*span_div)::setting_angle,setting_angle0 !取り付け角 [deg]
    real(8),dimension(0:2*span_div)::dihedral_angle,dihedral_angle0 !上反角 [deg]
    real(8),dimension(0:2*span_div)::foil1_mixture,foil2_mixture !翼型配合率
    real(8),dimension(0:2*span_div)::EIx,GJ !曲げ剛性，ねじれ定数
    real(8),dimension(0:2*span_div-1)::chord_cp !コントロールポイントにおけるコード長 [m]
    real(8),dimension(0:2*span_div)::deflection,theta !たわみ [m]，たわみ角 [deg]
    real(8),dimension(0:2*span_div)::phi !ねじれ角 [deg]
    real(8),dimension(0:2,0:2*span_div-1)::cp !コントロールポイント [m]
    real(8),dimension(0:2*span_div)::x,y,z !翼素座標 [m]
    real(8)::S,b,AR !翼面積 [m^2]，スパン [m]，アスペクト比 [-]
    !翼の空力計算
    real(8),dimension(0:14,0:1)::Cl_coef,Cdp_coef,Cm_coef !揚力係数，抗力係数，ピッチングモーメント係数を多項式近似したときの係数
    real(8)::dynamic_pressure !動圧 [N/m^2]
    real(8),dimension(0:2*span_div-1)::alpha_induced,alpha_effective !吹きおろし角 [deg]，有効迎角 [deg]
    real(8),dimension(0:2*span_div-1)::Re !Re数
    real(8),dimension(0:2*span_div-1)::Cl,Cdp,a0,a1 !局所揚力係数，有害抗力係数，N軸方向の係数，T軸方向の係数,断面揚力傾斜
    real(8),dimension(0:2*span_div-1)::Cm,Cm_cg !空力中心回り，重心位置まわりのピッチングモーメント係数
    real(8)::Cla !二次元揚力傾斜
    real(8),dimension(0:2*span_div-1)::Circulation,Circulation_old !循環Γ [m^2/s]
    real(8),dimension(0:2*span_div-1)::wi !吹きおろし速度 [m/s]
    real(8),dimension(0:2*span_div-1,0:2*span_div-1)::Qij
    !翼にはたらく力，モーメント
    real(8)::Induced_Drag !誘導抗力
    real(8),dimension(0:2*span_div-1)::dL,dDp,dW,dN,dT,dM,dM_cg !翼素揚力，有害抗力，重量，N軸方向の力，T軸方向の力,重心まわりのピッチングモーメント
    real(8),dimension(0:2*span_div)::Bending_Moment,Torque,Shear_Force !曲げモーメント [N*m]，ねじりモーメント [N*m]，せんだん力 [N]
    !地面効果の計算
    real(8),dimension(0:2*span_div-1,0:2*span_div-1)::yp,zp,ym,zm,ymp,zmp,Rpij,Rmij,Rpijm,Rmijm
    !収束判定用
    real(8)::Lift_wing,Lift_old,Induced_Drag_old
    !一般
    real(8)::sum,tmp,integral
    real(8)::pi,deg,rad !円周率，degへの変換，radへの変換
    integer::num,num1,num2,num3
    !カウンター
    integer::i,j,k,iteration
    !計算時間計測
    real(8)::CPU_time_s
    integer::initial_time,finish_time,time_rate,CPU_time_min !initial time, finish time, time rate

    call system_clock(initial_time) !開始時間の読み込み

    pi=4.0D0*atan(1.0D0) !円周率
    deg=(180.0D0/pi) !degへの変換
    rad=(pi/180.0D0) !radへの変換

    !値の読み込み
    !左翼端が0，右翼端が2*span_div
    open(10,file="Wing_data.txt") !ファイルを開く
        read(10,*) span_div !スパン方向分割数，コード方向分割数
        read(10,*) rho,mu !空気密度 [kg/m^3]，動粘性係数 [m^2/s]
        read(10,*) hac,hspar !空力中心位置
        read(10,*) dy !翼素間隔 [m]
        ds=dy*0.5D0
        read(10,*) chord_mac,hcg !平均空力翼弦長 [m]，重心位置 [-]
        read(10,*) zc(0:8)
        read(10,*) chord(0:span_div) !コード長 [m]
        read(10,*) Cl_coef(0:14,0) !揚力係数を多項式近似したときの係数
        read(10,*) Cdp_coef(0:14,0) !抗力係数を多項式近似したときの係数
        read(10,*) Cm_coef(0:14,0) !ピッチングモーメント係数を多項式近似したときの係数
        read(10,*) Cl_coef(0:14,1) !揚力係数を多項式近似したときの係数
        read(10,*) Cdp_coef(0:14,1) !抗力係数を多項式近似したときの係数
        read(10,*) Cm_coef(0:14,1) !ピッチングモーメント係数を多項式近似したときの係数
        read(10,*) setting_angle0(0:span_div) !取り付け角 [deg]
        read(10,*) dihedral_angle0(0:span_div) !上反角 [deg]
        read(10,*) foil1_mixture(0:span_div) !翼型配合率
        read(10,*) foil2_mixture(0:span_div) !翼型配合率
        read(10,*) EIx(0:span_div) !曲げ剛性
        read(10,*) GJ(0:span_div) !ねじれ定数
        read(10,*) dW(0:span_div) !翼素重量 [N]
    close(10) !ファイルを閉じる

    !左翼の値を右翼に格納
    !左翼端が0，右翼端が2*span_div
    do i=0,span_div
        chord(2*span_div-i)=chord(i)
        setting_angle0(2*span_div-i)=setting_angle0(i)
        dihedral_angle0(2*span_div-i)=dihedral_angle0(i)
        foil1_mixture(2*span_div-i)=foil1_mixture(i)
        foil2_mixture(2*span_div-i)=foil2_mixture(i)
        EIx(2*span_div-i)=EIx(i)
        GJ(2*span_div-i)=GJ(i)
        dW(2*span_div-i-1)=dW(i)
    end do

    !chord_cp,dS,Sを計算
    S=0.0D0
    do i=0,2*span_div-1
        chord_cp(i)=(chord(i)+chord(i+1))/2
        !dS(i)=chord_cp(i)*dy
        S=S+chord_cp(i)*dy
    end do

    !配列の初期化
    alpha_induced=0.0D0
    setting_angle=0.0D0
    dihedral_angle=0.0D0
    phi=0.0D0
    theta=0.0D0
    deflection=0.0D0
    circulation=0.0D0
    circulation_old=0.0D0

    !LLTの反復計算
    do iteration=0,iteration_max

        !左翼端が0，右翼端が2*span_div
        x=0.0D0;    y=0.0D0;    z=0.0D0
        do i=0,span_div !片翼分のサーフェスポイントのループ
            !上反角，取り付け角を更新
            !左翼
            setting_angle(i)=setting_angle0(i)+phi(i)
            dihedral_angle(i)=-dihedral_angle0(i)-theta(i) !左翼は負の値
            !右翼
            setting_angle(span_div+i)=setting_angle0(span_div-i)+phi(span_div-i)
            dihedral_angle(span_div+i)=dihedral_angle0(span_div-i)+theta(span_div-i)
            !y,zを計算
            if(i>0) then 
                !左翼
                y(span_div-i)=y(span_div-i+1)-dy*cos(-rad*dihedral_angle(span_div-i))
                z(span_div-i)=z(span_div-i+1)+dy*sin(-rad*dihedral_angle(span_div-i))
                !右翼
                y(span_div+i)=y(span_div+i-1)+dy*Cos(rad*dihedral_angle(span_div+i))
                z(span_div+i)=z(span_div+i-1)+dy*Sin(rad*dihedral_angle(span_div+i))
            end if
        end do

        !左翼端が0，右翼端が2*span_div
        cp=0.0D0
        do i=0,span_div-1 !片翼分のコントロールポイントのループ
            !cpの計算
            cp(0,i)=0.0D0
            cp(1,i)=(y(i)+y(i+1))/2.0D0
            cp(2,i)=(z(i)+z(i+1))/2.0D0
            num=2*span_div-1-i
            cp(0,num)=0.0D0
            cp(1,num)=(y(num)+y(num+1))/2.0D0
            cp(2,num)=(z(num)+z(num+1))/2.0D0
        end do

        !座標の計算
        yp=0.0D0;   zp=0.0D0;   ymp=0.0D0;  zmp=0.0D0
        do j=0,2*span_div-1
            do i=0,2*span_div-1
                !オリジナルの計算
                yp(i,j)=(cp(1,i)-cp(1,j))*Cos(rad*dihedral_angle(j)) &
                    +(cp(2,i)-cp(2,j))*Sin(rad*dihedral_angle(j))
                zp(i, j)=-(cp(1,i)-cp(1,j))*Sin(rad*dihedral_angle(j)) &
                    +(cp(2,i)-cp(2,j))*Cos(rad*dihedral_angle(j))
                !鏡像の計算
                ymp(i, j)=(cp(1,i)-cp(1,j))*Cos(-rad*dihedral_angle(j)) &
                    +(cp(2,i)-(-cp(2,j)-2.0D0*hE))*Sin(-rad*dihedral_angle(j))
                zmp(i,j)=-(cp(1,i)-cp(1,j))*Sin(-rad*dihedral_angle(j)) &
                    +(cp(2,i)-(-cp(2,j)-2.0D0*hE))*Cos(-rad*dihedral_angle(j))
        !    end do
        !end do
        
        !R+ij，R-ijの計算
        !do i=0,2*span_div-1
        !    do j=0,2*span_div-1
                !オリジナルの計算
                Rpij(i,j)=Sqrt((yp(i,j)-ds)**2+zp(i,j)**2)
                Rmij(i,j)=Sqrt((yp(i,j)+ds)**2+zp(i,j)**2)
                !鏡像の計算
                Rpijm(i,j)=Sqrt((ymp(i,j)-ds)**2+zmp(i,j)**2)
                Rmijm(i,j)=Sqrt((ymp(i,j)+ds)**2+zmp(i,j)**2)
        !    end do
        !end do
        
        !Qijの計算
        !do i=0,2*span_div-1
        !    do j=0,2*span_div-1
                Qij(i,j)=(-((yp(i,j)-ds)/(Rpij(i,j)*Rpij(i,j)))+((yp(i,j)+ds)/(Rmij(i,j)*Rmij(i,j)))) &
                            *Cos(rad*(dihedral_angle(i)-dihedral_angle(j))) &
                        +(-(zp(i,j)/(Rpij(i,j)*Rpij(i,j)))+(zp(i,j)/(Rmij(i,j)*Rmij(i,j)))) &
                            *Sin(rad*(dihedral_angle(i)-dihedral_angle(j))) &
                        +((ymp(i,j)-ds)/(Rpijm(i,j)*Rpijm(i,j))-(ymp(i,j)+ds)/(Rmijm(i,j)*Rmijm(i,j))) &
                            *Cos(rad*(dihedral_angle(i)+dihedral_angle(j))) &
                        +(zmp(i, j)/(Rpijm(i,j)*Rpijm(i,j))-zmp(i,j)/(Rmijm(i,j)*Rmijm(i,j))) &
                            *Sin(rad*(dihedral_angle(i)+dihedral_angle(j)))
            end do
        end do

        !循環Γの更新
        if(iteration>0) then 
            Circulation=Circulation_old+circulation_coef*(Circulation-Circulation_old)
        end if
        Circulation_old=Circulation
        !call write1D(2*span_div-1,circulation)

        !吹きおろしの計算
        !左翼端が0，右翼端が2*span_div
        epsilon=0.0D0
        wi=0.0D0
        do i=0,2*span_div-1
            !吹きおろし速度の計算
            do j=0,2*span_div-1
                wi(i)=wi(i)+(1.0D0/(4.0D0*pi))*Qij(i,j)*circulation(j)
            end do
            !誘導迎角の計算
            alpha_induced(i)=deg*atan(wi(i)/(Velocity-rad*r*cp(1,i)))
            epsilon=epsilon+alpha_induced(i)
        end do
        epsilon=epsilon/(2.0D0*real(span_div)-1.0D0) !主翼位置での平均吹きおろし角

        !左翼端が0，右翼端が2*span_div
        do i=0,span_div-1 !片翼分のコントロールポイントのループ
            !左翼の有効迎角の計算
            num=span_div-1-i
            alpha_effective(num)=alpha+setting_angle(num+1)-alpha_induced(num) & !有効迎角 [deg]＝全機迎角＋取り付け角ー誘導迎角
                +deg*atan((rad*p*cp(1,num))/(Velocity-rad*r*cp(1,num))) & !ロール角速度pによる影響を追加
                +deg*atan(Velocity*sin(rad*beta)*sin(rad*dihedral_angle(num+1))/(Velocity-rad*r*cp(1,num))) !横滑りによる影響を追加
            !右翼の有効迎角の計算
            num=span_div+i
            alpha_effective(num)=alpha+setting_angle(num)-alpha_induced(num) & !有効迎角 [deg]＝全機迎角＋取り付け角ー誘導迎角
                +deg*atan((rad*p*cp(1,num))/(Velocity-rad*r*cp(1,num))) & !ロール角速度による影響を追加
                +deg*atan(Velocity*sin(rad*beta)*sin(rad*dihedral_angle(num))/(Velocity-rad*r*cp(1,num))) !横滑りによる影響を追加
        end do
        do i=0,2*span_div-1
            !計算が発散しないようにalphaを制限する．
            if(alpha_effective(i)>alpha_max) alpha_effective(i)=alpha_max
            if(alpha_effective(i)<alpha_min) alpha_effective(i)=alpha_min
        end do

        !左翼端が0，右翼端が2*span_div
        do i=0,2*span_div-1 !片翼分のコントロールポイントのループ
            !レイノルズ数の計算
            Re(i)=((Velocity-rad*r*cp(1,i))*chord_cp(i))/mu !左翼
        end do
        do i=0,2*span_div-1
            !計算が発散しないようにRe数を制限する．
            if(Re(i)>Re_max) Re(i)=Re_max
            if(Re(i)<Re_min) Re(i)=Re_min
        end do

        !左翼端が0，右翼端が2*span_div
        Cl=0.0D0
        Cdp=0.0D0
        Cm=0.0D0
        !Drag=0.0D0
        !Force=0.0D0
        !Moment=0.0D0
        !Circulation_old=0.0D0
        do i=0,2*span_div-1 !両翼分のコントロールポイントのループ
            if(i<span_div) then
                num=i
            else
                num=i+1
            end if
            !断面揚力傾斜の計算 [1/deg]
            a0(i)=foil1_mixture(i)*Cl_coef(0,0)+foil2_mixture(i)*Cl_coef(0,1)
            a1(i)=foil1_mixture(i)*Cl_coef(1,0)+foil2_mixture(i)*Cl_coef(1,1)
            do j=0,8
                tmp=alpha_effective(i)**j
                Cl(i)=Cl(i)+foil1_mixture(num)*Cl_coef(j,0)*tmp+foil2_mixture(num)*Cl_coef(j,1)*tmp
                Cdp(i)=Cdp(i)+foil1_mixture(num)*Cdp_coef(j,0)*tmp+foil2_mixture(num)*Cdp_coef(j,1)*tmp
                Cm(i)=Cm(i)+foil1_mixture(num)*Cm_coef(j,0)*tmp+foil2_mixture(num)*Cm_coef(j,1)*tmp
            end do
            do j=9,14
                tmp=((Re(i)/100000.0D0)**(j-8))
                Cl(i)=Cl(i)+foil1_mixture(num)*Cl_coef(j,0)*tmp+foil2_mixture(num)*Cl_coef(j,1)*tmp
                Cdp(i)=Cdp(i)+foil1_mixture(num)*Cdp_coef(j,0)*tmp+foil2_mixture(num)*Cdp_coef(j,1)*tmp
                Cm(i)=Cm(i)+foil1_mixture(num)*Cm_coef(j,0)*tmp+foil2_mixture(num)*Cm_coef(j,1)*tmp
            end do
            Cm_cg(i)=Cm(i)+Cl(i)*(hcg-hac) !ピッチングモーメントを空力中心回りから桁位置まわりに変換

            !機体にはたらく空気力を計算
            tmp=chord_cp(i)*dy*Cos(rad*dihedral_angle(i)) !翼素面積
            dynamic_pressure=0.5D0*rho*(Velocity-rad*r*cp(1,i))**2 !動圧
            dL(i)=dynamic_pressure*tmp*CL(i) !翼素揚力
            dDp(i)=dynamic_pressure*tmp*Cdp(i) !翼素抗力
            dM_cg(i)=dynamic_pressure*tmp*chord_cp(i)*Cm_cg(i) !桁位置まわりの翼素ピッチングモーメント
            dM(i)=dynamic_pressure*tmp*chord_cp(i)*Cm(i) !空力中心まわりの翼素ピッチングモーメント
            
            !空気力を機体軸に対して分解
            !tmp=rad*(setting_angle(i)-alpha_effective(i))
            sum=rad*(alpha-alpha_induced(i)) &
                +atan((rad*p*cp(1,num))/(Velocity-rad*r*cp(1,num))) & !ロール角速度による影響を追加
                +atan(Velocity*sin(rad*beta)*sin(rad*dihedral_angle(num))/(Velocity-rad*r*cp(1,num))) !横滑りによる影響を追加
            !Cn(i)=CL(i)*cos(tmp)-Cdp(i)*sin(tmp)
            !Ct(i)=-CL(i)*sin(tmp)+Cdp(i)*cos(tmp)
            dN(i)=dynamic_pressure*tmp*(CL(i)*cos(sum)-Cdp(i)*sin(sum))
            dT(i)=dynamic_pressure*tmp*(-CL(i)*sin(sum)+Cdp(i)*cos(sum))
        end do

        !揚力，誘導抗力の計算
        Lift_wing=0.0D0
        Induced_Drag=0.0D0
        do i=0,2*span_div-1
            Lift_wing=Lift_wing+rho*(Velocity-rad*r*cp(1,i))*circulation(i)*dy*cos(rad*dihedral_angle(i))
            Induced_Drag=Induced_Drag+rho*wi(i)*circulation(i)*dy
        end do
        !write(*,*) Lift_wing

        !循環Γを計算
        !左翼端が0，右翼端が2*span_div
        do i=0,2*span_div-1 !両翼分のコントロールポイントのループ
            Circulation(i)=0.5D0*chord_cp(i)*(Velocity-rad*r*cp(1,i))*Cl(i)
        end do

        !曲げモーメント，ねじりモーメントを計算
        !左翼端が0，右翼端が2*span_div
        Bending_Moment=0.0D0
        Torque=0.0D0
        Shear_Force=0.0D0
        do i=1,span_div !翼端から翼根に向けてのループ
            dz=0.0D0
            do j=1,i !翼端からi番目の要素に向けてのループ

                !上反角がないとして曲げモーメントを計算．j>i．num1<num2
                num1=2*span_div-i
                num2=2*span_div-j
                !Bending_Moment(i)=Bending_Moment(i)+(dN(j-1)-dW(j-1))*(abs(y(j)-y(i))+0.5D0*dy) !左翼
                !Bending_Moment(num1)=Bending_Moment(num1)+(dN(num2)-dW(num2))*(y(num2)-y(num1)+0.5D0*dy) !右翼
                Bending_Moment(i)=Bending_Moment(i) &
                    +(dN(j-1)*cos(rad*dihedral_angle(j-1))-dW(j-1))*Abs(cp(1,j)-cp(1,i)) &
                    +dN(j-1)*sin(rad*dihedral_angle(j-1))*Abs(cp(2,j)-cp(2,j)) !左翼
                Bending_Moment(num1)=Bending_Moment(num1) &
                    +(dN(num2)*cos(rad*dihedral_angle(num2))-dW(num2))*Abs(cp(1,num2)-cp(1,num1)) &
                    +dN(num2)*sin(rad*dihedral_angle(num2))*Abs(cp(2,num2)-cp(2,num1)) !右翼

                !ねじりモーメント
                do k=i,j,-1 !i番目の要素からj番目の要素に向けてのループ．i>j
                    num3=2*span_div-k
                    dz(j)=dz(j)+dy*sin(rad*abs(dihedral_angle(k)-dihedral_angle(i))) !左翼
                    dz(num2)=dz(num2)+dy*sin(rad*(dihedral_angle(num3)-dihedral_angle(num1))) !右翼
                end do
                Torque(i)=Torque(i)+dM_cg(j-1)+dz(j)*dT(j-1) !左翼
                Torque(num1)=Torque(num1)+dM_cg(num2)+dz(num2)*dT(num2) !右翼

                !剪断力
                Shear_Force(i)=Shear_Force(i)+(dN(j-1)-dW(j-1))
                Shear_Force(num1)=Shear_Force(num1)+(dN(num2)-dW(num2))

            end do
        end do
        Bending_Moment(span_div)=0.5D0*Bending_Moment(span_div)
        Torque(span_div)=0.50D0*Torque(span_div)
        Shear_Force(span_div)=0.5D0*Shear_Force(span_div)
        
        !たわみ，たわみ角，ねじれ角の計算
        !左翼端が0，右翼端が2*span_div
        phi=0.0D0
        theta=0.0D0
        deflection=0.0D0
        do i=span_div-1,0,-1 !翼中心から翼端へのループ
            num=2*span_div-1-i
            !write(*,*) num+1
            !たわみ角
            theta(i)=theta(i+1)+Bending_Moment(i)*dy/EIx(i) !左翼
            theta(num+1)=theta(num)+Bending_Moment(num)*dy/EIx(num+1)  !右翼
            !たわみ
            deflection(i)=deflection(i+1)+theta(i+1)*dy !左翼
            deflection(num+1)=deflection(num)+theta(num)*dy !右翼
            !ねじれ角
            phi(i)=phi(i+1)+(dM(i)*(dy/2.0D0))/GJ(i+1)+(Torque(i+1)*dy)/GJ(i+1) !左翼
            phi(num+1)=phi(num)+(dM(num)*(dy/2.0D0))/GJ(num)+(Torque(num)*dy)/GJ(num) !右翼
        end do
        theta=deg*theta ![deg]への変換
        phi=deg*phi ![deg]への変換

        if(iteration>1) then
            if(error>abs((Lift_wing-Lift_old)/Lift_old)) then
                exit
            else
                Lift_old=Lift_wing
            end if
        end if
        !write(*,*) Lift_wing

        call system_clock(finish_time,time_rate) !終了時間の読み込み
        CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
        !if(CPU_time_s>1.000D0) exit !解析時間が1秒を超えたらおそらく発散しているので強制終了．

    end do

    Lift=Lift+Lift_wing
    Drag=Drag+Induced_Drag

    do i=0,2*span_div-1 !両翼分のコントロールポイントのループ
        !空気力X,Y,Z,モーメントL,M,Nを計算
        Drag=Drag+dDp(i) !有害抗力
        Force(0)=Force(0)-dT(i) !X [N]
        IF(i<span_div) then
            Force(1)=Force(1)+dN(i)*sin(rad*dihedral_angle(i+1)) !Y [N]
            Force(2)=Force(2)-dN(i)*cos(rad*dihedral_angle(i+1)) !Z [N]
            Moment(0)=Moment(0)-dN(i)*cos(rad*dihedral_angle(i+1))*cp(1,i) !L [N*m]
        else
            Force(1)=Force(1)-dN(i)*sin(rad*dihedral_angle(i)) !Y [N]
            Force(2)=Force(2)-dN(i)*cos(rad*dihedral_angle(i)) !Z [N]
            Moment(0)=Moment(0)-dN(i)*cos(rad*dihedral_angle(i))*cp(1,i) !L [N*m]
        end if
        !Force(2)=Force(2)-dN(i) !Z [N]
        !Moment(0)=Moment(0)-dN(i)*cp(1,i) !L [N*m]
        Moment(1)=Moment(1)+dM(i)+dT(i)*cp(2,i) !M[N*m] 空力中心周り
        Moment(2)=Moment(2)+dT(i)*cp(1,i) !N [N*m]
    end do
    Moment(1)=Moment(1)+Lift_wing*chord_mac*(hcg-dh-hac) !重心位置周り

end subroutine

subroutine Wing_calculation_linear(Lift,Drag,Force,Moment,Velocity,alpha,beta,p,r,hE,dh,epsilon,span_div)
    !機体の姿勢角による翼素の高度への影響は無視．横滑りによって渦糸が斜めに出る影響は無視．
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan,sqrt,sin,cos,tan !プログラム中で使用する関数の名前を宣言
    !引数
    real(8),intent(INOUT)::Lift !揚力 [N]
    real(8),intent(INOUT)::Drag !抗力 [N]
    real(8),dimension(0:2),intent(INOUT)::Force !機体にはたらく外力 [X,Y,Z] [N]
    real(8),dimension(0:2),intent(INOUT)::Moment !機体にはたらくモーメント [L,M,N] [N*m]
    real(8),intent(IN)::Velocity,alpha,beta,p,r,hE,dh !対気速度 [m/s]，迎角 [deg]，横滑り角 [deg]，角速度 [deg/s]高度 [m]，重心移動距離 [-]
    real(8),intent(OUT)::epsilon !主翼位置での吹きおろし角 [deg]
    integer,intent(INOUT)::span_div
    !LLT
    !パラメータ
    integer::iteration_max=1 !最大反復計算回数
    real(8),parameter::alpha_max=15.0D0 !計算を発散させないためのalphaの最大値
    real(8),parameter::alpha_min=-5.0D0 !計算を発散させないためのalphaの最小値
    real(8),parameter::Re_max =1000000.0D0 !計算を発散させないためのRe数の最小値
    real(8),parameter::Re_min =100000.0D0 !計算を発散させないためのRe数の最小値
    !環境変数
    real(8)::rho,mu !空気密度 [kg/m^3]，動粘性係数 [m^2/s]
    !解析用
    real(8)::chord_div !スパン方向分割数，コード方向分割数
    !翼の幾何学形状
    real(8)::dy,ds !翼素間隔，パネルの反復 [m]
    real(8),dimension(0:2*span_div-1)::dz !高さの差 [m]
    real(8)::hac,hspar !空力中心位置，桁位置
    real(8)::hcg !重心位置
    real(8)::chord_mac !平均空力翼弦長 [m]
    integer::num_mac !平均空力翼弦長になるリブ番号
    real(8),dimension(0:8)::zc !中心線を多項式近似したときの係数
    real(8),dimension(0:2*span_div)::chord !コード長 [m]
    real(8),dimension(0:2*span_div)::setting_angle,setting_angle0 !取り付け角 [deg]
    real(8),dimension(0:2*span_div)::dihedral_angle,dihedral_angle0 !上反角 [deg]
    real(8),dimension(0:2*span_div)::foil1_mixture,foil2_mixture !翼型配合率
    real(8),dimension(0:2*span_div)::EIx,GJ !曲げ剛性，ねじれ定数
    real(8),dimension(0:2*span_div-1)::chord_cp !コントロールポイントにおけるコード長 [m]
    real(8),dimension(0:2*span_div)::deflection,theta !たわみ [m]，たわみ角 [deg]
    real(8),dimension(0:2*span_div)::phi !ねじれ角 [deg]
    real(8),dimension(0:2,0:2*span_div-1)::cp !コントロールポイント [m]
    real(8),dimension(0:2*span_div)::x,y,z !翼素座標 [m]
    real(8)::S,b,AR !翼面積 [m^2]，スパン [m]，アスペクト比 [-]
    !翼の空力計算
    real(8),dimension(0:14,0:1)::Cl_coef,Cdp_coef,Cm_coef !揚力係数，抗力係数，ピッチングモーメント係数を多項式近似したときの係数
    real(8)::dynamic_pressure !動圧 [N/m^2]
    real(8),dimension(0:2*span_div-1)::alpha_induced,alpha_effective !吹きおろし角 [deg]，有効迎角 [deg]
    real(8),dimension(0:2*span_div-1)::Re !Re数
    real(8),dimension(0:2*span_div-1)::Cl,Cdp,a0,a1 !局所揚力係数，有害抗力係数，N軸方向の係数，T軸方向の係数,断面揚力傾斜
    real(8),dimension(0:2*span_div-1)::Cm,Cm_cg !空力中心回り，重心位置まわりのピッチングモーメント係数
    real(8)::Cla !二次元揚力傾斜
    real(8),dimension(0:2*span_div-1)::Circulation,Circulation_old !循環Γ [m^2/s]
    real(8),dimension(0:2*span_div-1)::wi !吹きおろし速度 [m/s]
    real(8),dimension(0:2*span_div-1,0:2*span_div-1)::Qij
    !翼にはたらく力，モーメント
    real(8)::Lift_wing,Induced_Drag !誘導抗力
    real(8),dimension(0:2*span_div-1)::dL,dDp,dW,dN,dT,dM,dM_cg !翼素揚力，有害抗力，重量，N軸方向の力，T軸方向の力,重心まわりのピッチングモーメント
    real(8),dimension(0:2*span_div)::Bending_Moment,Torque,Shear_Force !曲げモーメント [N*m]，ねじりモーメント [N*m]，せんだん力 [N]
    !地面効果の計算
    real(8),dimension(0:2*span_div-1,0:2*span_div-1)::yp,zp,ym,zm,ymp,zmp,Rpij,Rmij,Rpijm,Rmijm
    !連立方程式用
    real(8),dimension(0:2*span_div-1,0:2*span_div-1)::matrix1
    real(8),dimension(0:2*span_div-1)::matrix2,answer
    !一般
    real(8)::sum,tmp,integral
    real(8)::pi,deg,rad !円周率，degへの変換，radへの変換
    integer::num,num1,num2,num3
    !カウンター
    integer::i,j,k,iteration
    !計算時間計測
    real(8)::CPU_time_s
    integer::initial_time,finish_time,time_rate,CPU_time_min !initial time, finish time, time rate

    call system_clock(initial_time) !開始時間の読み込み

    pi=4.0D0*atan(1.0D0) !円周率
    deg=(180.0D0/pi) !degへの変換
    rad=(pi/180.0D0) !radへの変換

    !値の読み込み
    !左翼端が0，右翼端が2*span_div
    open(10,file="Wing_data.txt") !ファイルを開く
        read(10,*) span_div !スパン方向分割数，コード方向分割数
        read(10,*) rho,mu !空気密度 [kg/m^3]，動粘性係数 [m^2/s]
        read(10,*) hac,hspar !空力中心位置
        read(10,*) dy !翼素間隔 [m]
        ds=dy*0.5D0
        read(10,*) chord_mac,hcg !平均空力翼弦長 [m]，重心位置 [-]
        read(10,*) zc(0:8)
        read(10,*) chord(0:span_div) !コード長 [m]
        read(10,*) Cl_coef(0:14,0) !揚力係数を多項式近似したときの係数
        read(10,*) Cdp_coef(0:14,0) !抗力係数を多項式近似したときの係数
        read(10,*) Cm_coef(0:14,0) !ピッチングモーメント係数を多項式近似したときの係数
        read(10,*) Cl_coef(0:14,1) !揚力係数を多項式近似したときの係数
        read(10,*) Cdp_coef(0:14,1) !抗力係数を多項式近似したときの係数
        read(10,*) Cm_coef(0:14,1) !ピッチングモーメント係数を多項式近似したときの係数
        read(10,*) setting_angle0(0:span_div) !取り付け角 [deg]
        read(10,*) dihedral_angle0(0:span_div) !上反角 [deg]
        read(10,*) foil1_mixture(0:span_div) !翼型配合率
        read(10,*) foil2_mixture(0:span_div) !翼型配合率
        read(10,*) EIx(0:span_div) !曲げ剛性
        read(10,*) GJ(0:span_div) !ねじれ定数
        read(10,*) dW(0:span_div) !翼素重量 [N]
    close(10) !ファイルを閉じる

    !左翼の値を右翼に格納
    !左翼端が0，右翼端が2*span_div
    do i=0,span_div
        chord(2*span_div-i)=chord(i)
        setting_angle0(2*span_div-i)=setting_angle0(i)
        dihedral_angle0(2*span_div-i)=dihedral_angle0(i)
        foil1_mixture(2*span_div-i)=foil1_mixture(i)
        foil2_mixture(2*span_div-i)=foil2_mixture(i)
        EIx(2*span_div-i)=EIx(i)
        GJ(2*span_div-i)=GJ(i)
        dW(2*span_div-i-1)=dW(i)
    end do

    !chord_cp,dS,Sを計算
    S=0.0D0
    do i=0,2*span_div-1
        chord_cp(i)=(chord(i)+chord(i+1))/2
        !dS(i)=chord_cp(i)*dy
        S=S+chord_cp(i)*dy
    end do

    !配列の初期化
    alpha_induced=0.0D0
    setting_angle=0.0D0
    dihedral_angle=0.0D0
    phi=0.0D0
    theta=0.0D0
    deflection=0.0D0
    circulation=0.0D0
    circulation_old=0.0D0

    !LLTの反復計算
    do iteration=0,iteration_max

        !左翼端が0，右翼端が2*span_div
        x=0.0D0;    y=0.0D0;    z=0.0D0
        do i=0,span_div !片翼分のサーフェスポイントのループ
            !上反角，取り付け角を更新
            !左翼
            setting_angle(i)=setting_angle0(i)+phi(i)
            dihedral_angle(i)=-dihedral_angle0(i)-theta(i) !左翼は負の値
            !右翼
            setting_angle(span_div+i)=setting_angle0(span_div-i)+phi(span_div-i)
            dihedral_angle(span_div+i)=dihedral_angle0(span_div-i)+theta(span_div-i)
            !y,zを計算
            if(i>0) then 
                !左翼
                y(span_div-i)=y(span_div-i+1)-dy*cos(-rad*dihedral_angle(span_div-i))
                z(span_div-i)=z(span_div-i+1)+dy*sin(-rad*dihedral_angle(span_div-i))
                !右翼
                y(span_div+i)=y(span_div+i-1)+dy*Cos(rad*dihedral_angle(span_div+i))
                z(span_div+i)=z(span_div+i-1)+dy*Sin(rad*dihedral_angle(span_div+i))
            end if
        end do

        !左翼端が0，右翼端が2*span_div
        cp=0.0D0
        do i=0,span_div-1 !片翼分のコントロールポイントのループ
            !cpの計算
            cp(0,i)=0.0D0
            cp(1,i)=(y(i)+y(i+1))/2.0D0
            cp(2,i)=(z(i)+z(i+1))/2.0D0
            num=2*span_div-1-i
            cp(0,num)=0.0D0
            cp(1,num)=(y(num)+y(num+1))/2.0D0
            cp(2,num)=(z(num)+z(num+1))/2.0D0
        end do

        !座標の計算
        yp=0.0D0;   zp=0.0D0;   ymp=0.0D0;  zmp=0.0D0
        do j=0,2*span_div-1
            do i=0,2*span_div-1
                !オリジナルの計算
                yp(i,j)=(cp(1,i)-cp(1,j))*Cos(rad*dihedral_angle(j)) &
                    +(cp(2,i)-cp(2,j))*Sin(rad*dihedral_angle(j))
                zp(i, j)=-(cp(1,i)-cp(1,j))*Sin(rad*dihedral_angle(j)) &
                    +(cp(2,i)-cp(2,j))*Cos(rad*dihedral_angle(j))
                !鏡像の計算
                ymp(i, j)=(cp(1,i)-cp(1,j))*Cos(-rad*dihedral_angle(j)) &
                    +(cp(2,i)-(-cp(2,j)-2.0D0*hE))*Sin(-rad*dihedral_angle(j))
                zmp(i,j)=-(cp(1,i)-cp(1,j))*Sin(-rad*dihedral_angle(j)) &
                    +(cp(2,i)-(-cp(2,j)-2.0D0*hE))*Cos(-rad*dihedral_angle(j))
        !    end do
        !end do
        
        !R+ij，R-ijの計算
        !do i=0,2*span_div-1
        !    do j=0,2*span_div-1
                !オリジナルの計算
                Rpij(i,j)=Sqrt((yp(i,j)-ds)**2+zp(i,j)**2)
                Rmij(i,j)=Sqrt((yp(i,j)+ds)**2+zp(i,j)**2)
                !鏡像の計算
                Rpijm(i,j)=Sqrt((ymp(i,j)-ds)**2+zmp(i,j)**2)
                Rmijm(i,j)=Sqrt((ymp(i,j)+ds)**2+zmp(i,j)**2)
        !    end do
        !end do
        
        !Qijの計算
        !do i=0,2*span_div-1
        !    do j=0,2*span_div-1
                Qij(i,j)=(-((yp(i,j)-ds)/(Rpij(i,j)*Rpij(i,j)))+((yp(i,j)+ds)/(Rmij(i,j)*Rmij(i,j)))) &
                            *Cos(rad*(dihedral_angle(i)-dihedral_angle(j))) &
                        +(-(zp(i,j)/(Rpij(i,j)*Rpij(i,j)))+(zp(i,j)/(Rmij(i,j)*Rmij(i,j)))) &
                            *Sin(rad*(dihedral_angle(i)-dihedral_angle(j))) &
                        +((ymp(i,j)-ds)/(Rpijm(i,j)*Rpijm(i,j))-(ymp(i,j)+ds)/(Rmijm(i,j)*Rmijm(i,j))) &
                            *Cos(rad*(dihedral_angle(i)+dihedral_angle(j))) &
                        +(zmp(i, j)/(Rpijm(i,j)*Rpijm(i,j))-zmp(i,j)/(Rmijm(i,j)*Rmijm(i,j))) &
                            *Sin(rad*(dihedral_angle(i)+dihedral_angle(j)))
            end do
        end do

        !左翼端が0，右翼端が2*span_div
        do i=0,span_div-1 !片翼分のコントロールポイントのループ
            !左翼の有効迎角の計算
            num=span_div-1-i
            alpha_effective(num)=alpha+setting_angle(num+1) & !有効迎角 [deg]＝全機迎角＋取り付け角
                +deg*atan((rad*p*cp(1,num))/(Velocity-rad*r*cp(1,num))) & !ロール角速度pによる影響を追加
                +deg*atan(Velocity*sin(rad*beta)*sin(rad*dihedral_angle(num+1))/(Velocity-rad*r*cp(1,num))) !横滑りによる影響を追加
            !右翼の有効迎角の計算
            num=span_div+i
            alpha_effective(num)=alpha+setting_angle(num) & !有効迎角 [deg]＝全機迎角＋取り付け角
                +deg*atan((rad*p*cp(1,num))/(Velocity-rad*r*cp(1,num))) & !ロール角速度による影響を追加
                +deg*atan(Velocity*sin(rad*beta)*sin(rad*dihedral_angle(num))/(Velocity-rad*r*cp(1,num))) !横滑りによる影響を追加
        end do
        do i=0,2*span_div-1
            !計算が発散しないようにalphaを制限する．
            if(alpha_effective(i)>alpha_max) alpha_effective(i)=alpha_max
            if(alpha_effective(i)<alpha_min) alpha_effective(i)=alpha_min
        end do

        !左翼端が0，右翼端が2*span_div
        Cl=0.0D0
        do i=0,2*span_div-1 !両翼分のコントロールポイントのループ
            if(i<span_div) then
                num=i
            else
                num=i+1
            end if
            !断面揚力傾斜の計算 [1/deg]
            a0(i)=foil1_mixture(i)*Cl_coef(0,0)+foil2_mixture(i)*Cl_coef(0,1)
            a1(i)=foil1_mixture(i)*Cl_coef(1,0)+foil2_mixture(i)*Cl_coef(1,1)
            do j=0,8
                tmp=alpha_effective(i)**j
                Cl(i)=Cl(i)+foil1_mixture(num)*Cl_coef(j,0)*tmp+foil2_mixture(num)*Cl_coef(j,1)*tmp
            end do
        end do

        !行列A，Bを計算
        do j=0,2*span_div-1
            do i=0,2*span_div-1
                !行列Aの計算
                If(i==j) Then
                    matrix1(i,j)=(8.0D0*pi)/(chord_cp(i)*a1(i))+Qij(i,j)
                Else
                    matrix1(i,j)=0.0D0+Qij(i,j)
                End If
            end do
            !行列Bの計算
            matrix2(j)=4.0D0*pi*(Velocity-rad*r*cp(1,j))*CL(j)/a1(j)
        end do
        
        !境界条件の方程式を解いて循環を求める
        call Ans(2*span_div,matrix1,matrix2,answer)
        do i=0,2*span_div-1
            circulation(i)=answer(i)
        end do
        
        !吹きおろしの計算
        !左翼端が0，右翼端が2*span_div
        epsilon=0.0D0
        wi=0.0D0
        !Cl=0.0D0
        !Cdp=0.0D0
        !Cm=0.0D0
        do i=0,2*span_div-1
            !吹きおろし速度の計算
            do j=0,2*span_div-1
                wi(i)=wi(i)+(1.0D0/(4.0D0*pi))*Qij(i,j)*circulation(j)
            end do
            !誘導迎角の計算
            alpha_induced(i)=deg*atan(wi(i)/(Velocity-rad*r*cp(1,i)))
            epsilon=epsilon+alpha_induced(i)
        end do
        epsilon=epsilon/(2.0D0*real(span_div)-1.0D0) !主翼位置での平均吹きおろし角

        !局所揚力係数CLから有効迎角αeを求める
        do i=0,2*span_div-1
            alpha_effective(i)=alpha_effective(i)-alpha_induced(i)
        end do

        !左翼端が0，右翼端が2*span_div
        do i=0,2*span_div-1 !片翼分のコントロールポイントのループ
            !レイノルズ数の計算
            Re(i)=((Velocity-rad*r*cp(1,i))*chord_cp(i))/mu !左翼
        end do
        do i=0,2*span_div-1
            !計算が発散しないようにRe数を制限する．
            if(Re(i)>Re_max) Re(i)=Re_max
            if(Re(i)<Re_min) Re(i)=Re_min
        end do

        !左翼端が0，右翼端が2*span_div
        Cl=0.0D0
        Cdp=0.0D0
        Cm=0.0D0
        do i=0,2*span_div-1 !両翼分のコントロールポイントのループ
            if(i<span_div) then
                num=i
            else
                num=i+1
            end if
            !断面揚力傾斜の計算 [1/deg]
            a0(i)=foil1_mixture(i)*Cl_coef(0,0)+foil2_mixture(i)*Cl_coef(0,1)
            a1(i)=foil1_mixture(i)*Cl_coef(1,0)+foil2_mixture(i)*Cl_coef(1,1)
            do j=0,8
                tmp=alpha_effective(i)**j
                Cl(i)=Cl(i)+foil1_mixture(num)*Cl_coef(j,0)*tmp+foil2_mixture(num)*Cl_coef(j,1)*tmp
                Cdp(i)=Cdp(i)+foil1_mixture(num)*Cdp_coef(j,0)*tmp+foil2_mixture(num)*Cdp_coef(j,1)*tmp
                Cm(i)=Cm(i)+foil1_mixture(num)*Cm_coef(j,0)*tmp+foil2_mixture(num)*Cm_coef(j,1)*tmp
            end do
            do j=9,14
                tmp=((Re(i)/100000.0D0)**(j-8))
                Cl(i)=Cl(i)+foil1_mixture(num)*Cl_coef(j,0)*tmp+foil2_mixture(num)*Cl_coef(j,1)*tmp
                Cdp(i)=Cdp(i)+foil1_mixture(num)*Cdp_coef(j,0)*tmp+foil2_mixture(num)*Cdp_coef(j,1)*tmp
                Cm(i)=Cm(i)+foil1_mixture(num)*Cm_coef(j,0)*tmp+foil2_mixture(num)*Cm_coef(j,1)*tmp
            end do
            Cm_cg(i)=Cm(i)+Cl(i)*(hcg-hac) !ピッチングモーメントを空力中心回りから桁位置まわりに変換

            !機体にはたらく空気力を計算
            tmp=chord_cp(i)*dy*Cos(rad*dihedral_angle(i)) !翼素面積
            dynamic_pressure=0.5D0*rho*(Velocity-rad*r*cp(1,i))**2 !動圧
            dL(i)=dynamic_pressure*tmp*CL(i) !翼素揚力
            dDp(i)=dynamic_pressure*tmp*Cdp(i) !翼素抗力
            dM_cg(i)=dynamic_pressure*tmp*chord_cp(i)*Cm_cg(i) !桁位置まわりの翼素ピッチングモーメント
            dM(i)=dynamic_pressure*tmp*chord_cp(i)*Cm(i) !空力中心まわりの翼素ピッチングモーメント
            
            !空気力を機体軸に対して分解
            !tmp=rad*(setting_angle(i)-alpha_effective(i))
            sum=rad*(alpha-alpha_induced(i)) &
                +atan((rad*p*cp(1,num))/(Velocity-rad*r*cp(1,num))) & !ロール角速度による影響を追加
                +atan(Velocity*sin(rad*beta)*sin(rad*dihedral_angle(num))/(Velocity-rad*r*cp(1,num))) !横滑りによる影響を追加
            !Cn(i)=CL(i)*cos(tmp)-Cdp(i)*sin(tmp)
            !Ct(i)=-CL(i)*sin(tmp)+Cdp(i)*cos(tmp)
            dN(i)=dynamic_pressure*tmp*(CL(i)*cos(sum)-Cdp(i)*sin(sum))
            dT(i)=dynamic_pressure*tmp*(-CL(i)*sin(sum)+Cdp(i)*cos(sum))
        end do

        !循環Γを計算
        !左翼端が0，右翼端が2*span_div
        do i=0,2*span_div-1 !両翼分のコントロールポイントのループ
            Circulation(i)=0.5D0*chord_cp(i)*(Velocity-rad*r*cp(1,i))*Cl(i)
        end do
        !epsilon=epsilon/(2.0D0*real(span_div)-1.0D0) !主翼位置での平均吹きおろし角

        !揚力，誘導抗力の計算
        Lift_wing=0.0D0
        Induced_Drag=0.0D0
        do i=0,2*span_div-1
            Lift_wing=Lift_wing+rho*(Velocity-rad*r*cp(1,i))*circulation(i)*dy*cos(rad*dihedral_angle(i))
            Induced_Drag=Induced_Drag+rho*wi(i)*circulation(i)*dy
        end do
        !write(*,*) Lift_wing

        !曲げモーメント，ねじりモーメントを計算
        !左翼端が0，右翼端が2*span_div
        Bending_Moment=0.0D0
        Torque=0.0D0
        Shear_Force=0.0D0
        do i=1,span_div !翼端から翼根に向けてのループ
            dz=0.0D0
            do j=1,i !翼端からi番目の要素に向けてのループ

                !上反角がないとして曲げモーメントを計算．j>i．num1<num2
                num1=2*span_div-i
                num2=2*span_div-j
                !Bending_Moment(i)=Bending_Moment(i)+(dN(j-1)-dW(j-1))*(abs(y(j)-y(i))+0.5D0*dy) !左翼
                !Bending_Moment(num1)=Bending_Moment(num1)+(dN(num2)-dW(num2))*(y(num2)-y(num1)+0.5D0*dy) !右翼
                Bending_Moment(i)=Bending_Moment(i) &
                    +(dN(j-1)*cos(rad*dihedral_angle(j-1))-dW(j-1))*Abs(cp(1,j)-cp(1,i)) &
                    +dN(j-1)*sin(rad*dihedral_angle(j-1))*Abs(cp(2,j)-cp(2,j)) !左翼
                Bending_Moment(num1)=Bending_Moment(num1) &
                    +(dN(num2)*cos(rad*dihedral_angle(num2))-dW(num2))*Abs(cp(1,num2)-cp(1,num1)) &
                    +dN(num2)*sin(rad*dihedral_angle(num2))*Abs(cp(2,num2)-cp(2,num1)) !右翼

                !ねじりモーメント
                do k=i,j,-1 !i番目の要素からj番目の要素に向けてのループ．i>j
                    num3=2*span_div-k
                    dz(j)=dz(j)+dy*sin(rad*abs(dihedral_angle(k)-dihedral_angle(i))) !左翼
                    dz(num2)=dz(num2)+dy*sin(rad*(dihedral_angle(num3)-dihedral_angle(num1))) !右翼
                end do
                Torque(i)=Torque(i)+dM_cg(j-1)+dz(j)*dT(j-1) !左翼
                Torque(num1)=Torque(num1)+dM_cg(num2)+dz(num2)*dT(num2) !右翼

                !剪断力
                Shear_Force(i)=Shear_Force(i)+(dN(j-1)-dW(j-1))
                Shear_Force(num1)=Shear_Force(num1)+(dN(num2)-dW(num2))

            end do
        end do
        Bending_Moment(span_div)=0.5D0*Bending_Moment(span_div)
        Torque(span_div)=0.50D0*Torque(span_div)
        Shear_Force(span_div)=0.5D0*Shear_Force(span_div)
        
        !たわみ，たわみ角，ねじれ角の計算
        !左翼端が0，右翼端が2*span_div
        phi=0.0D0
        theta=0.0D0
        deflection=0.0D0
        do i=span_div-1,0,-1 !翼中心から翼端へのループ
            num=2*span_div-1-i
            !write(*,*) num+1
            !たわみ角
            theta(i)=theta(i+1)+Bending_Moment(i)*dy/EIx(i) !左翼
            theta(num+1)=theta(num)+Bending_Moment(num)*dy/EIx(num+1)  !右翼
            !たわみ
            deflection(i)=deflection(i+1)+theta(i+1)*dy !左翼
            deflection(num+1)=deflection(num)+theta(num)*dy !右翼
            !ねじれ角
            phi(i)=phi(i+1)+(dM(i)*(dy/2.0D0))/GJ(i+1)+(Torque(i+1)*dy)/GJ(i+1) !左翼
            phi(num+1)=phi(num)+(dM(num)*(dy/2.0D0))/GJ(num)+(Torque(num)*dy)/GJ(num) !右翼
        end do
        theta=deg*theta ![deg]への変換
        phi=deg*phi ![deg]への変換

        call system_clock(finish_time,time_rate) !終了時間の読み込み
        CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
        !if(CPU_time_s>1.000D0) exit !解析時間が1秒を超えたらおそらく発散しているので強制終了．

    end do

    Lift=Lift+Lift_wing
    Drag=Drag+Induced_Drag

    do i=0,2*span_div-1 !両翼分のコントロールポイントのループ
        !空気力X,Y,Z,モーメントL,M,Nを計算
        Drag=Drag+dDp(i) !有害抗力
        Force(0)=Force(0)-dT(i) !X [N]
        IF(i<span_div) then
            Force(1)=Force(1)+dN(i)*sin(rad*dihedral_angle(i+1)) !Y [N]
            Force(2)=Force(2)-dN(i)*cos(rad*dihedral_angle(i+1)) !Z [N]
            Moment(0)=Moment(0)-dN(i)*cos(rad*dihedral_angle(i+1))*cp(1,i) !L [N*m]
        else
            Force(1)=Force(1)-dN(i)*sin(rad*dihedral_angle(i)) !Y [N]
            Force(2)=Force(2)-dN(i)*cos(rad*dihedral_angle(i)) !Z [N]
            Moment(0)=Moment(0)-dN(i)*cos(rad*dihedral_angle(i))*cp(1,i) !L [N*m]
        end if
        !Force(2)=Force(2)-dN(i) !Z [N]
        !Moment(0)=Moment(0)-dN(i)*cp(1,i) !L [N*m]
        Moment(1)=Moment(1)+dM(i)+dT(i)*cp(2,i) !M[N*m] 空力中心周り
        Moment(2)=Moment(2)+dT(i)*cp(1,i) !N [N*m]
    end do
    Moment(1)=Moment(1)+Lift_wing*chord_mac*(hcg-dh-hac) !重心位置周り
    
end subroutine

subroutine Tail_calculation(Lift,Drag,Force,Moment,Velocity,alpha,q,de,dh,epsilon,span_div,chord_div)
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan,sqrt,sin,cos,tan !プログラム中で使用する関数の名前を宣言
    !引数
    real(8),intent(INOUT)::Lift !揚力 [N]
    real(8),intent(INOUT)::Drag !抗力 [N]
    real(8),dimension(0:2),intent(INOUT)::Force !機体にはたらく外力 [X,Y,Z] [N]
    real(8),dimension(0:2),intent(INOUT)::Moment !機体にはたらくモーメント [L,M,N] [N*m]
    real(8),intent(IN)::Velocity,alpha,q,de,dh !対気速度 [m/s]，迎角 [deg]，角速度 [deg/s]，エレベータ舵角 [deg]，重心移動距離 [-]
    real(8),intent(IN)::epsilon !主翼位置での吹きおろし角 [deg]
    integer,intent(INOUT)::span_div,chord_div
    !VLM
    !パラメータ
    real(8),parameter::alpha_max=10.0D0 !計算を発散させないためのalphaの最大値
    real(8),parameter::alpha_min=-10.0D0 !計算を発散させないためのalphaの最小値
    real(8),parameter::Re_max =1000000.0D0 !計算を発散させないためのRe数の最小値
    real(8),parameter::Re_min =100000.0D0 !計算を発散させないためのRe数の最小値
    !環境変数
    real(8)::rho,mu !空気密度 [kg/m^3]，動粘性係数 [m^2/s]
    !翼の幾何学形状
    real(8)::dy !翼素間隔 [m]
    real(8)::hac,hspar !空力中心位置，桁位置
    real(8)::hcg !重心位置
    real(8)::chord_mac !平均空力翼弦長 [m]
    real(8)::lt !重心-水平尾翼中心間距離 [m]
    real(8)::it !水平尾翼取り付け角 [deg]
    real(8),dimension(0:8)::zc !中心線を多項式近似したときの係数
    real(8),dimension(0:2*span_div)::chord !コード長 [m]
    real(8),dimension(0:2*span_div-1)::chord_cp !コントロールポイントにおけるコード長 [m]
    real(8),dimension(0:2*span_div-1)::dS !翼素面積 [m^2]
    real(8),dimension(0:2,0:2*span_div*chord_div-1)::cp !コントロールポイント [m]
    real(8),dimension(0:2,0:2*span_div*chord_div-1)::sp1,sp2 !サーフェスポイント [m]
    real(8),dimension(0:(2*span_div+1)*chord_div+2*span_div)::x,y,z !翼素座標 [m]
    real(8),dimension(0:2*span_div*chord_div-1)::dzdx !キャンバーラインの傾き
    real(8),dimension(0:2*span_div*chord_div-1,0:2*span_div*chord_div-1)::Cij !j番目のU字渦がi番目のcpに誘導する速度を表す係数
    real(8)::dx1,dx2,dy1,dy2 !dx1=xi-x1j,dx2=xi-x2j,dy1=yi-y1j,dy2=yi-y2j
    real(8)::r1,r2 !r1=sqrt(dx1^2+dy1^2),r2=sqrt(dx2^2+dy2^2)
    real(8)::S,b,AR !翼面積 [m^2]，スパン [m]，アスペクト比 [-]
    real(8)::x0,y0,z0
    !翼の空力計算
    real(8)::alpha_tail !尾翼迎角
    real(8),dimension(0:14)::Cl_coef,Cdp_coef,Cm_coef !揚力係数，抗力係数，ピッチングモーメント係数を多項式近似したときの係数
    real(8)::dynamic_pressure !動圧 [N/m^2]
    real(8),dimension(0:2*span_div-1)::alpha_induced,alpha_effective !吹きおろし角 [deg]，有効迎角 [deg]
    real(8),dimension(0:2*span_div-1)::Re !Re数
    real(8),dimension(0:2*span_div-1)::Cl,Cdp,Cn,Ct !局所揚力係数，有害抗力係数，N軸方向の係数，T軸方向の係数
    real(8),dimension(0:2*span_div-1)::Cm,Cm_cg !空力中心回り，重心位置まわりのピッチングモーメント係数
    real(8)::Cla !二次元揚力傾斜
    real(8),dimension(0:2*span_div*chord_div-1)::Circulation !循環Γ [m^2/s]
    real(8),dimension(0:2*span_div-1)::Circulation_chord !循環Γ [m^2/s]
    real(8),dimension(0:2*span_div-1)::wi !吹きおろし速度 [m/s]
    real(8)::wi_x,wi_y,wi_z !各方向の吹きおろし速度 [m/s]
    !翼にはたらく力，モーメント
    real(8)::Lift_tail !尾翼揚力
    real(8)::Induced_Drag !誘導抗力
    real(8),dimension(0:2*span_div-1)::dL,dDp,dW,dN,dT,dM !翼素揚力，有害抗力，重量，N軸方向の力，T軸方向の力,重心まわりのピッチングモーメント
    !一般
    real(8)::sum,tmp,integral
    real(8)::pi,deg,rad !円周率，degへの変換，radへの変換
    real(8),dimension(:,:),allocatable::matrix1
    real(8),dimension(:),allocatable::matrix2
    integer::num,num1,num2,num3
    !カウンター
    integer(8)::i,j,k,iteration

    pi=4.0D0*atan(1.0D0) !円周率
    deg=(180.0D0/pi) !degへの変換
    rad=(pi/180.0D0) !radへの変換

    !値の読み込み
    !左翼端が0，右翼端が2*span_div
    open(10,file="Tail_data.txt") !ファイルを開く
        read(10,*) span_div,chord_div !スパン方向分割数，コード方向分割数
        read(10,*) rho,mu !空気密度 [kg/m^3]，動粘性係数 [m^2/s]
        read(10,*) hac,hspar !空力中心位置
        read(10,*) dy !翼素間隔 [m]
        read(10,*) chord_mac,lt,it !平均空力翼弦長 [m]，重心-水平尾翼中心間距離 [m]，水平尾翼取り付け角 [deg]
        read(10,*) zc(0:8)
        read(10,*) chord(0:span_div) !コード長 [m]
        read(10,*) Cl_coef(0:14) !揚力係数を多項式近似したときの係数
        read(10,*) Cdp_coef(0:14) !抗力係数を多項式近似したときの係数
        read(10,*) Cm_coef(0:14) !ピッチングモーメント係数を多項式近似したときの係数
    close(10) !ファイルを閉じる

    !尾翼迎角を計算
    alpha_tail=alpha+it+de-2*epsilon+deg*atan(rad*q*(lt+chord_mac*dh)/Velocity)

    !左翼の値を右翼に格納
    !左翼端が0，右翼端が2*span_div
    do i=0,span_div
        chord(2*span_div-i)=chord(i)
    end do

    !chord_cp,dS,Sを計算
    S=0.0D0
    do i=0,2*span_div-1
        chord_cp(i)=(chord(i)+chord(i+1))/2
        dS(i)=chord_cp(i)*dy
        S=S+dS(i)
    end do

    !x,y,zを計算
    !桁位置が原点
    !左翼端が0，右翼端が2*span_div
    z=0.0D0
    do i=0,2*span_div
        do j=0,chord_div
            num=i*chord_div+i+j
            If(i<span_div) Then
                x(num)=-chord(i)*hspar+chord(i)*(real(j)/chord_div)
                y(num)=-dy*(real(span_div)-real(i))
                do k=0,8
                    z(num)=z(num)+chord(i)*zc(k)*((real(j)/chord_div)**k)
                end do
            Else
                x(num)=-chord(2*span_div-i)*hspar+chord(2*span_div-i)*(real(j)/chord_div)
                y(num)=dy*(real(i)-real(span_div))
                do k=0,8
                    z(num)=z(num)+chord(2*span_div-i)*zc(k)*((real(j)/chord_div)**k)
                end do
            End If
        end do
    end do

    !sp,cp,dz/dxを計算
    do i = 0,2*span_div-1
        do j =0,chord_div-1
            num1=i*chord_div+j !cp
            num2=i*chord_div+i+j !xyz
            !平面翼なのでz=0
            cp(0,num1)=((x(num2)*1.0D0+x(num2+1)*3.0D0)/4.0D0+(x(num2+chord_div+1)*1.0D0+x(num2+chord_div+2)*3.0D0)/4.0D0)/2.0D0
            cp(1,num1)=((y(num2)*1.0D0+y(num2+1)*3.0D0)/4.0D0+(y(num2+chord_div+1)*1.0D0+y(num2+chord_div+2)*3.0D0)/4.0D0)/2.0D0
            cp(2,num1)=0.0D0
            sp1(0,num1)=(x(num2)*3.0D0+x(num2+1)*1.0D0)/4.0D0
            sp1(1,num1)=(y(num2)*3.0D0+y(num2+1)*1.0D0)/4.0D0
            sp1(2,num1)=0.0D0
            sp2(0,num1)=(x(num2+chord_div+1)*3.0D0+x(num2+chord_div+2)*1.0D0)/4.0D0
            sp2(1,num1)=(y(num2+chord_div+1)*3.0D0+y(num2+chord_div+2)*1.0D0)/4.0D0
            sp2(2,num1)=0.0D0
            dzdx(num1)=((z(num2+1)+z(num2+chord_div+2))/2.0D0-(z(num2)+z(num2+chord_div+1))/2.0D0) / &
                ((x(num2+1)+x(num2+chord_div+2))/2.0D0-(x(num2)+x(num2+chord_div+1))/2.0D0)
            !write(*,'(I4)') num1
        end do
    end do

    !Cijを計算
    do i=0,2*span_div*chord_div-1
        do j=0,2*span_div*chord_div-1
            dx1=cp(0,i)-sp1(0,j) !dx1=xi-x1j
            dx2=cp(0,i)-sp2(0,j) !dx2=xi-x2j
            dy1=cp(1,i)-sp1(1,j) !dy1=yi-y1j
            dy2=cp(1,i)-sp2(1,j) !dy2=yi-y2j
            r1=Sqrt(dx1**2+dy1**2) !r1=sqrt(dx1^2+dy1^2)
            r2=Sqrt(dx2**2+dy2**2) !r2=sqrt(dx2^2+dy2^2)
            Cij(i,j)=(1.0D0/(dx1*dy2-dy1*dx2))*(((dx1-dx2)*dx1+(dy1-dy2)*dy1)/r1 &
                -((dx1-dx2)*dx2+(dy1-dy2)*dy2)/r2)-(1.0D0/dy1)*(1.0D0+dx1/r1)+(1.0D0/dy2)*(1.0D0+dx2/r2)
        end do
    end do

    !境界条件の方程式を計算する
    allocate(matrix1(0:2*span_div*chord_div-1,0:2*span_div*chord_div-1)) 
    allocate(matrix2(0:2*span_div*chord_div-1))
    do i=0,2*span_div*chord_div-1
        do j=0,2*span_div*chord_div-1
            matrix1(i, j) = Cij(i, j)
        end do
    end do
    do i=0,2*span_div*chord_div-1
        matrix2(i)=-4.0D0*pi*Velocity*Sin(rad*alpha_tail-dzdx(i))
    end do

    !境界条件の方程式を解いて循環を求める
    circulation=0.0D0
    call Ans(2*span_div*chord_div,matrix1,matrix2,circulation)

    !循環のコード方向の和から，揚力係数，局所有効迎角，誘導迎角，吹きおろし速度を計算する
    circulation_chord=0.0D0
    do i=0,2*span_div-1
        do j=0,chord_div-1
            num1=i*chord_div+j !cp
            circulation_chord(i)=circulation_chord(i)+ circulation(num1)
        end do
    end do
        
    !循環分布から吹きおろしを計算する
    wi=0.0D0
    do i=0,2*span_div-1
        do j=0,2*span_div-1 !左翼の循環による吹きおろし
            wi(i)=wi(i)+(1.0D0/(4.0D0*pi)) &
                *(1.0D0/(sp1(1,j*chord_div)-cp(1,i*chord_div))-1.0D0/(sp2(1,j*chord_div)-cp(1,i*chord_div)))*circulation_chord(j)
        end do
        alpha_induced(i)=deg*Atan(wi(i)/Velocity)
        alpha_effective(i)=alpha_tail+alpha_induced(i)
        Re(i)=(Velocity*chord_cp(i))/mu !左翼
    end do
    do i=0,2*span_div-1
        !計算が発散しないようにRe数を制限する．
        if(Re(i)>Re_max) Re(i)=Re_max
        if(Re(i)<Re_min) Re(i)=Re_min
        !計算が発散しないようにalphaを制限する．
        if(alpha_effective(i)>alpha_max) alpha_effective(i)=alpha_max
        if(alpha_effective(i)<alpha_min) alpha_effective(i)=alpha_min
    end do

    !揚力，誘導抗力を計算する
    !wは下向きが正
    Induced_Drag=0.0D0
    Lift_tail=0.0D0
    do i=0,2*span_div-1
        Lift_tail=Lift_tail+rho*Velocity*circulation_chord(i)*dy
        Induced_Drag=Induced_Drag-rho*wi(i)*circulation_chord(i)*dy
    end do
    Lift=Lift+Lift_tail

    !左翼端が0，右翼端が2*span_div
    Cl=0.0D0
    Cdp=0.0D0
    Cm=0.0D0
    do i=0,2*span_div-1 !両翼分のコントロールポイントのループ
        if(i<span_div) then
            num=i
        else
            num=i+1
        end if
        do j=0,8
            tmp=alpha_effective(i)**j
            Cl(i)=Cl(i)+Cl_coef(j)*tmp
            Cdp(i)=Cdp(i)+Cdp_coef(j)*tmp
            Cm(i)=Cm(i)+Cm_coef(j)*tmp
        end do
        do j=9,14
            tmp=((Re(i)/100000.0D0)**(j-8))
            Cl(i)=Cl(i)+Cl_coef(j)*tmp
            Cdp(i)=Cdp(i)+Cdp_coef(j)*tmp
            Cm(i)=Cm(i)+Cm_coef(j)*tmp
        end do
        Cm_cg(i)=Cm(i)+Cl(i)*(hspar-hac) !ピッチングモーメントを空力中心回りから桁位置まわりに変換

        !機体にはたらく空気力を計算
        dynamic_pressure=0.5D0*rho*Velocity**2 !動圧
        dL(i)=dynamic_pressure*dS(i)*CL(i) !翼素揚力
        dDp(i)=dynamic_pressure*dS(i)*Cdp(i) !翼素抗力
        dM(i)=dynamic_pressure*dS(i)*chord_cp(i)*Cm_cg(i) !桁位置まわりの翼素ピッチングモーメント
        
        !空気力を機体軸に対して分解
        tmp=rad*(alpha_tail-alpha_induced(i))
        Cn(i)=CL(i)*cos(tmp)-Cdp(i)*sin(tmp)
        Ct(i)=-CL(i)*sin(tmp)+Cdp(i)*cos(tmp)
        dN(i)=dynamic_pressure*dS(i)*Cn(i)
        dT(i)=dynamic_pressure*dS(i)*Ct(i)

        !空気力X,Y,Z,モーメントL,M,Nを計算
        Drag=Drag+dDp(i) !有害抗力
        Force(0)=Force(0)-dT(i) !X [N]
        Force(2)=Force(2)-dN(i) !Z [N]
        Moment(1)=Moment(1)+dM(i) !M [N*m]

    end do
    Moment(1)=Moment(1)-Lift_tail*(lt+chord_mac*dh)

end subroutine

subroutine Fin_calculation(Lift,Drag,Force,Moment,Velocity,beta,p,r,dr,dh,span_div,chord_div)
    !胴体の影響は無視
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan,sqrt,sin,cos,tan !プログラム中で使用する関数の名前を宣言
    !引数
    real(8),intent(INOUT)::Lift !揚力 [N]
    real(8),intent(INOUT)::Drag !抗力 [N]
    real(8),dimension(0:2),intent(INOUT)::Force !機体にはたらく外力 [X,Y,Z] [N]
    real(8),dimension(0:2),intent(INOUT)::Moment !機体にはたらくモーメント [L,M,N] [N*m]
    real(8),intent(IN)::Velocity,beta,p,r,dr,dh !対気速度 [m/s]，横滑り角 [deg]，角速度 [deg/s]，舵角 [deg]，重心移動距離 [-]
    integer,intent(INOUT)::span_div,chord_div
    !VLM
    !パラメータ
    real(8),parameter::Fin_efficient=0.800D0 !垂直尾翼効率
    real(8),parameter::alpha_max=10.0D0 !計算を発散させないためのalphaの最大値
    real(8),parameter::alpha_min=-10.0D0 !計算を発散させないためのalphaの最小値
    real(8),parameter::Re_max =1000000.0D0 !計算を発散させないためのRe数の最小値
    real(8),parameter::Re_min =100000.0D0 !計算を発散させないためのRe数の最小値
    !環境変数
    real(8)::rho,mu !空気密度 [kg/m^3]，動粘性係数 [m^2/s]
    !翼の幾何学形状
    real(8)::dy !翼素間隔 [m]
    real(8)::hac,hspar !空力中心位置，桁位置
    real(8)::hcg !重心位置
    real(8)::chord_mac !平均空力翼弦長 [m]
    real(8)::lf !重心-垂直尾翼中心間距離 [m]
    real(8)::zf !重心-垂直尾翼下端距離 [m]
    real(8)::yf !垂直尾翼面積重心位置 [m]
    real(8),dimension(0:8)::zc !中心線を多項式近似したときの係数
    real(8),dimension(0:span_div)::chord !コード長 [m]
    real(8),dimension(0:span_div-1)::chord_cp !コントロールポイントにおけるコード長 [m]
    real(8),dimension(0:span_div-1)::dS !翼素面積 [m^2]
    real(8),dimension(0:2,0:span_div*chord_div-1)::cp !コントロールポイント [m]
    real(8),dimension(0:2,0:span_div*chord_div-1)::sp1,sp2 !サーフェスポイント [m]
    real(8),dimension(0:(span_div+1)*chord_div+span_div)::x,y,z !翼素座標 [m]
    real(8),dimension(0:span_div*chord_div-1)::dzdx !キャンバーラインの傾き
    real(8),dimension(0:span_div*chord_div-1,0:span_div*chord_div-1)::Cij !j番目のU字渦がi番目のcpに誘導する速度を表す係数
    real(8)::dx1,dx2,dy1,dy2 !dx1=xi-x1j,dx2=xi-x2j,dy1=yi-y1j,dy2=yi-y2j
    real(8)::r1,r2 !r1=sqrt(dx1^2+dy1^2),r2=sqrt(dx2^2+dy2^2)
    real(8)::S,b,AR !翼面積 [m^2]，スパン [m]，アスペクト比 [-]
    real(8)::x0,y0,z0
    !翼の空力計算
    real(8)::alpha_fin !尾翼迎角
    real(8),dimension(0:14)::Cl_coef,Cdp_coef,Cm_coef !揚力係数，抗力係数，ピッチングモーメント係数を多項式近似したときの係数
    real(8)::dynamic_pressure !動圧 [N/m^2]
    real(8),dimension(0:span_div-1)::alpha_induced,alpha_effective !吹きおろし角 [deg]，有効迎角 [deg]
    real(8),dimension(0:span_div-1)::Re !Re数
    real(8),dimension(0:span_div-1)::Cl,Cdp,Cn,Ct !局所揚力係数，有害抗力係数，N軸方向の係数，T軸方向の係数
    real(8),dimension(0:span_div-1)::Cm,Cm_cg !空力中心回り，重心位置まわりのピッチングモーメント係数
    real(8),dimension(0:span_div*chord_div-1)::Circulation !循環Γ [m^2/s]
    real(8),dimension(0:span_div-1)::Circulation_chord !循環Γ [m^2/s]
    real(8),dimension(0:span_div-1)::wi !吹きおろし速度 [m/s]
    real(8)::wi_x,wi_y,wi_z !各方向の吹きおろし速度 [m/s]
    !翼にはたらく力，モーメント
    real(8)::Lift_Fin !尾翼揚力
    real(8)::Induced_Drag !誘導抗力
    real(8),dimension(0:span_div-1)::dL,dDp,dW,dN,dT,dM !翼素揚力，有害抗力，重量，N軸方向の力，T軸方向の力,重心まわりのピッチングモーメント
    !一般
    real(8)::sum,tmp,integral
    real(8)::pi,deg,rad !円周率，degへの変換，radへの変換
    real(8),dimension(:,:),allocatable::matrix1
    real(8),dimension(:),allocatable::matrix2
    integer::num,num1,num2,num3
    !カウンター
    integer(8)::i,j,k,iteration



    pi=4.0D0*atan(1.0D0) !円周率
    deg=(180.0D0/pi) !degへの変換
    rad=(pi/180.0D0) !radへの変換

    !値の読み込み
    !下端が0，上端がspan_div
    open(10,file="Fin_data.txt") !ファイルを開く
        read(10,*) span_div,chord_div !スパン方向分割数，コード方向分割数
        read(10,*) rho,mu !空気密度 [kg/m^3]，動粘性係数 [m^2/s]
        read(10,*) hac,hspar !空力中心位置
        read(10,*) dy !翼素間隔 [m]
        read(10,*) chord_mac,lf,zf !平均空力翼弦長 [m]，重心-水平尾翼中心間距離 [m]，重心-垂直尾翼下端距離 [m]
        read(10,*) zc(0:8)
        read(10,*) chord(0:span_div) !コード長 [m]
        read(10,*) Cl_coef(0:14) !揚力係数を多項式近似したときの係数
        read(10,*) Cdp_coef(0:14) !抗力係数を多項式近似したときの係数
        read(10,*) Cm_coef(0:14) !ピッチングモーメント係数を多項式近似したときの係数
    close(10) !ファイルを閉じる

    !垂直尾翼迎角を計算
    alpha_fin=beta-dr-deg*atan(rad*r*(lf+chord_mac*dh)/Velocity)

    !chord_cp,dS,Sを計算
    S=0.0D0
    do i=0,span_div-1
        chord_cp(i)=(chord(i)+chord(i+1))/2
        dS(i)=chord_cp(i)*dy
        S=S+dS(i)
    end do

    !x,y,zを計算
    !桁位置が原点
    !下端が0，上端がspan_div
    z=0.0D0
    do i=0,span_div
        do j=0,chord_div
            num=i*chord_div+i+j
            x(num)=-chord(i)*hspar+chord(i)*(real(j)/chord_div)
            y(num)=dy*real(i)
            do k=0,8
                z(num)=z(num)+chord(i)*zc(k)*((real(j)/chord_div)**k)
            end do
            if(isnan(x(i))) then
                write(*,'(A12,I12,I12)') "x",iteration,i
                stop
            end if
            if(isnan(y(i))) then
                write(*,'(A12,I12,I12)') "y",iteration,i
                stop
            end if
            if(isnan(z(i))) then
                write(*,'(A12,I12,I12)') "z",iteration,i
                stop
            end if
        end do
    end do

    !sp,cp,dz/dxを計算
    do i = 0,span_div-1
        do j =0,chord_div-1
            num1=i*chord_div+j !cp
            num2=i*chord_div+i+j !xyz
            !平面翼なのでz=0
            cp(0,num1)=((x(num2)*1.0D0+x(num2+1)*3.0D0)/4.0D0+(x(num2+chord_div+1)*1.0D0+x(num2+chord_div+2)*3.0D0)/4.0D0)/2.0D0
            cp(1,num1)=((y(num2)*1.0D0+y(num2+1)*3.0D0)/4.0D0+(y(num2+chord_div+1)*1.0D0+y(num2+chord_div+2)*3.0D0)/4.0D0)/2.0D0
            cp(2,num1)=0.0D0
            sp1(0,num1)=(x(num2)*3.0D0+x(num2+1)*1.0D0)/4.0D0
            sp1(1,num1)=(y(num2)*3.0D0+y(num2+1)*1.0D0)/4.0D0
            sp1(2,num1)=0.0D0
            sp2(0,num1)=(x(num2+chord_div+1)*3.0D0+x(num2+chord_div+2)*1.0D0)/4.0D0
            sp2(1,num1)=(y(num2+chord_div+1)*3.0D0+y(num2+chord_div+2)*1.0D0)/4.0D0
            sp2(2,num1)=0.0D0
            dzdx(num1)=((z(num2+1)+z(num2+chord_div+2))/2.0D0-(z(num2)+z(num2+chord_div+1))/2.0D0) / &
                ((x(num2+1)+x(num2+chord_div+2))/2.0D0-(x(num2)+x(num2+chord_div+1))/2.0D0)
            if(isnan(dzdx(num1))) then
                write(*,'(A12,I12,I12)') "dzdx",num1
                stop
            end if
        end do
    end do

    !Cijを計算
    do i=0,span_div*chord_div-1
        do j=0,span_div*chord_div-1
            dx1=cp(0,i)-sp1(0,j) !dx1=xi-x1j
            dx2=cp(0,i)-sp2(0,j) !dx2=xi-x2j
            dy1=cp(1,i)-sp1(1,j) !dy1=yi-y1j
            dy2=cp(1,i)-sp2(1,j) !dy2=yi-y2j
            r1=Sqrt(dx1**2+dy1**2) !r1=sqrt(dx1^2+dy1^2)
            r2=Sqrt(dx2**2+dy2**2) !r2=sqrt(dx2^2+dy2^2)
            Cij(i,j)=(1.0D0/(dx1*dy2-dy1*dx2))*(((dx1-dx2)*dx1+(dy1-dy2)*dy1)/r1 &
                -((dx1-dx2)*dx2+(dy1-dy2)*dy2)/r2)-(1.0D0/dy1)*(1.0D0+dx1/r1)+(1.0D0/dy2)*(1.0D0+dx2/r2)
        end do
    end do

    !境界条件の方程式を計算する
    allocate(matrix1(0:span_div*chord_div-1,0:span_div*chord_div-1)) 
    allocate(matrix2(0:span_div*chord_div-1))
    do i=0,span_div*chord_div-1
        do j=0,span_div*chord_div-1
            matrix1(i, j) = Cij(i, j)
        end do
    end do
    do i=0,span_div*chord_div-1
        matrix2(i)=-4.0D0*pi*Velocity*Sin(rad*alpha_fin-dzdx(i))
        if(isnan(matrix2(i))) then
            write(*,'(A12,I12,I12)') "matrix2",i
            stop
        end if
    end do

    !境界条件の方程式を解いて循環を求める
    circulation=0.0D0
    call Ans(span_div*chord_div,matrix1,matrix2,circulation)

    !循環のコード方向の和から，揚力係数，局所有効迎角，誘導迎角，吹きおろし速度を計算する
    circulation_chord=0.0D0
    do i=0,span_div-1
        do j=0,chord_div-1
            num1=i*chord_div+j !cp
            circulation_chord(i)=circulation_chord(i)+ circulation(num1)
            if(isnan(circulation_chord(i))) then
                write(*,'(A12,I12,I12)') "circulation_chord",i
                stop
            end if
        end do
    end do
        
    !循環分布から吹きおろしを計算する
    wi=0.0D0
    do i=0,span_div-1
        do j=0,span_div-1 !左翼の循環による吹きおろし
            wi(i)=wi(i)+(1.0D0/(4.0D0*pi)) &
                *(1.0D0/(sp1(1,j*chord_div)-cp(1,i*chord_div))-1.0D0/(sp2(1,j*chord_div)-cp(1,i*chord_div)))*circulation_chord(j)
        end do
        alpha_induced(i)=deg*Atan(wi(i)/Velocity)
        alpha_effective(i)=alpha_fin+alpha_induced(i)
        Re(i)=(Velocity*chord_cp(i))/mu !左翼
        if(isnan(alpha_effective(i))) then
            write(*,'(A12,I12,I12)') "alpha_effective",i
            stop
        end if
        if(isnan(Re(i))) then
            write(*,'(A12,I12,I12)') "Re",i
            stop
        end if
    end do
    do i=0,span_div-1
        !計算が発散しないようにRe数を制限する．
        if(Re(i)>Re_max) Re(i)=Re_max
        if(Re(i)<Re_min) Re(i)=Re_min
        !計算が発散しないようにalphaを制限する．
        if(alpha_effective(i)>alpha_max) alpha_effective(i)=alpha_max
        if(alpha_effective(i)<alpha_min) alpha_effective(i)=alpha_min
    end do

    !揚力，誘導抗力を計算する
    !wは下向きが正
    Induced_Drag=0.0D0
    Lift_fin=0.0D0
    do i=0,span_div-1
        Lift_Fin=Lift_Fin+rho*Velocity*circulation_chord(i)*dy
        Induced_Drag=Induced_Drag-rho*wi(i)*circulation_chord(i)*dy
    end do
    Lift=Lift+Lift_Fin

    !片翼面積重心yf [m]の計算
    yf=0.0D0
    do i=0,span_div-1
        yf=yf+chord_cp(i)*cp(1,i*chord_div)*dy
    end do
    yf=(1.0D0/S)*yf

    !下端が0，上端がspan_div
    Cl=0.0D0
    Cdp=0.0D0
    Cm=0.0D0
    do i=0,span_div-1 !両翼分のコントロールポイントのループ
        num=i
        do j=0,8
            tmp=alpha_effective(i)**j
            Cl(i)=Cl(i)+Cl_coef(j)*tmp
            Cdp(i)=Cdp(i)+Cdp_coef(j)*tmp
            Cm(i)=Cm(i)+Cm_coef(j)*tmp
        end do
        do j=9,14
            tmp=((Re(i)/100000.0D0)**(j-8))
            Cl(i)=Cl(i)+Cl_coef(j)*tmp
            Cdp(i)=Cdp(i)+Cdp_coef(j)*tmp
            Cm(i)=Cm(i)+Cm_coef(j)*tmp
        end do
        Cm_cg(i)=Cm(i)+Cl(i)*(hspar-hac) !ピッチングモーメントを空力中心回りから桁位置まわりに変換

        !機体にはたらく空気力を計算
        dynamic_pressure=0.5D0*rho*Velocity**2 !動圧
        dL(i)=dynamic_pressure*dS(i)*CL(i) !翼素揚力
        dDp(i)=dynamic_pressure*dS(i)*Cdp(i) !翼素抗力
        dM(i)=dynamic_pressure*dS(i)*chord_cp(i)*Cm_cg(i) !桁位置まわりの翼素ピッチングモーメント
        
        !空気力を機体軸に対して分解
        tmp=rad*(alpha_fin-alpha_induced(i))
        Cn(i)=CL(i)*cos(tmp)-Cdp(i)*sin(tmp)
        Ct(i)=-CL(i)*sin(tmp)+Cdp(i)*cos(tmp)
        dN(i)=dynamic_pressure*dS(i)*Cn(i)
        dT(i)=dynamic_pressure*dS(i)*Ct(i)

        !空気力X,Y,Z,モーメントL,M,Nを計算
        Drag=Drag+dDp(i) !有害抗力
        Force(1)=Force(1)-dN(i) !Y [N]
        !Moment(2)=Moment(2)+dM(i) !N [N*m]

    end do
    Moment(0)=Moment(0)-(zf+yf)*Lift_fin*Fin_efficient !L [N*m]
    Moment(2)=Moment(2)+(lf+chord_mac*dh)*Lift_fin*Fin_efficient !N [N*m]

end subroutine

subroutine Rotate_bector(axis,theta,matrix1,matrix2)
    !3行1列のベクトルを指定された軸まわりに指定された角度だけ回転させる．
    implicit none
    integer,intent(IN)::axis !0:x軸，1：y軸，2：z軸
    real(8),intent(IN)::theta
    real(8),dimension(0:2,0:2)::rtt
    real(8),dimension(0:2,0:0),intent(IN)::matrix1
    real(8),dimension(0:2,0:0),intent(OUT)::matrix2
    real(8),dimension(0:2,0:0)::matrix
    real(8)::l1, m1, l2, m2, l3, m3
    !回転行列の定義
    l1 = cos(theta)
    m1 = sin(theta)

    if(axis==0) then !x軸まわりの回転行列
        rtt = 0.0D0
        rtt(0,0) = 1.0D0
        rtt(1,1) = l1
        rtt(2,1) = -m1
        rtt(1,2) = m1
        rtt(2,2) = l1
    elseif(axis==1) then !y軸まわりの回転行列
        rtt = 0.0D0
        rtt(0,0) = l1
        rtt(2,0) = m1
        rtt(1,1) = 1.0D0
        rtt(0,2) = -m1
        rtt(2,2) = l1
    elseif(axis==2) then !z軸まわりの回転行列
        rtt = 0.0D0
        rtt(0,0) = l1
        rtt(1,0) = -m1
        rtt(0,1) = m1
        rtt(1,1) = l1
        rtt(2,2) = 1.0D0
    end if

    matrix = matrix1
    call matrix_multiplication(3,3,rtt,3,1,matrix,matrix2) !軸まわりに回転

end subroutine

subroutine Transform_coordinate(phi,theta,psi,matrix1,matrix2)
    !3行1列のベクトルを絶対座標系から機体軸座標に変換する
    implicit none
    real(8),intent(IN)::phi,theta,psi
    real(8),dimension(0:2,0:2)::rtt
    real(8),dimension(0:2,0:0),intent(IN)::matrix1
    real(8),dimension(0:2,0:0),intent(OUT)::matrix2
    real(8),dimension(0:2,0:0)::matrix

    matrix=matrix1
    call rotate_bector(2,psi,matrix,matrix2) !ヨー方向に回転
    matrix=matrix2
    call rotate_bector(1,theta,matrix,matrix2) !ピッチ方向に回転
    matrix=matrix2
    call rotate_bector(0,phi,matrix,matrix2) !ロール方向に回転

end subroutine

subroutine Inverse_Transform_coordinate(phi,theta,psi,matrix1,matrix2)
    !3行1列のベクトルを機体座標系から絶対座標系に変換する
    implicit none
    real(8),intent(IN)::phi,theta,psi
    real(8),dimension(0:2,0:2)::rtt
    real(8),dimension(0:2,0:0),intent(IN)::matrix1
    real(8),dimension(0:2,0:0),intent(OUT)::matrix2
    real(8),dimension(0:2,0:0)::matrix

    matrix=matrix1
    call rotate_bector(0,-phi,matrix,matrix2) !ロール方向に逆回転
    matrix=matrix2
    call rotate_bector(1,-theta,matrix,matrix2) !ピッチ方向に逆回転
    matrix=matrix2
    call rotate_bector(2,-psi,matrix,matrix2) !ヨー方向に逆回転

end subroutine

subroutine matrix_multiplication(m,n,matrix1,p,q,matrix2,matrix3)
    !行列の積の計算
    implicit none
    integer,intent(IN)::m,n,p,q
    real(8),dimension(0:m-1,0:n-1),intent(IN)::matrix1
    real(8),dimension(0:p-1,0:q-1),intent(IN)::matrix2
    real(8),dimension(0:m-1,0:q-1),intent(OUT)::matrix3
    real(8)::a
    integer::i,j,k
    matrix3 = 0.0D0

    if(n /= p) then !行と列の数の不一致を知らせる
        write(*,*) "IMPOSSIBLE!"
    else !行列の積
        do i = 0, m-1
            do j = 0, q-1
                a = 0.0D0
                do k = 0, n-1
                a = a+matrix1(i,k)*matrix2(k,j)
                end do
                matrix3(i,j) = a
            end do
        end do
    end if

end subroutine matrix_multiplication

subroutine Ans(n,A,b,x)
    !ガウスの消去法によって，連立方程式Ax=bを解く．
    !Aはn×n行列，x，bはn×1行列
    !参考文献：「解析塾秘伝！有限要素法のつくり方」　p.82-
    !引数
    integer,intent(IN)::n
    real(8),dimension(0:n-1,0:n-1),intent(INOUT)::A
    real(8),dimension(0:n-1),intent(INOUT)::b
    real(8),dimension(0:n-1),intent(OUT)::x(0:n-1) !解を入れる配列
    !変数の宣言
    real(8)::pivot !マトリックスの体格成分
    real(8)::p !計算に使用するマトリックスの成分
    !カウンター
    integer::r,c,rr,cc

    !前進消去
    do r = 0,n - 1
        !対角成分をpivotに代入
        pivot = A(r, r)
        do c = r,N - 1    
            A(r, c) = A(r, c) / pivot
        end do
        b(r) = b(r) / pivot
        do rr = r + 1,n - 1
            p = A(rr, r)
            do cc = r,N - 1
                A(rr, cc) = A(rr, cc) - p * A(r, cc)
            end do
            b(rr) = b(rr) - p * b(r)
        end do
    end do

    !後退代入
    do r = N - 1,0,-1
        x(r) = b(r)
        do c = r + 1,N - 1
            x(r) = x(r) - A(r, c) * x(c)
        end do
    end do

End subroutine

subroutine write1D(n,x) !一次元配列の表示
    implicit none
    integer,intent(IN)::n
    real(8),dimension(0:n-1),intent(IN)::x
    integer::i

    do i = 0, n-1
    write(*,'("x(",I4,")=",F12.3)') i,x(i)
    end do

end subroutine

subroutine write2D(m,n,x) !二次元配列の表示
    implicit none
    integer,intent(IN)::m,n
    real(8),dimension(0:m-1,0:n-1),intent(IN)::x
    integer::i,j

    do i = 0, m-1
    do j = 0, n-1
        write(*,'("x(",I4,",",I2,")=",F12.6,4X)',advance='NO') i,j,x(i,j)
    end do
    write(*,*)
    end do

end subroutine
