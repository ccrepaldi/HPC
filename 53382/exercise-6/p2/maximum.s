	.cstring
LC2:
	.ascii "maximum.f90\0"
	.text
_MAIN__:
LFB0:
	pushq	%rbp
LCFI0:
	movq	%rsp, %rbp
LCFI1:
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$528, %rsp
LCFI2:
	movl	$1, %ebx
	movl	$20000000, %r12d
L3:
	testl	%r12d, %r12d
	jle	L2
	movslq	%ebx, %rax
	leaq	-1(%rax), %r14
	leal	-20000000(%rbx), %ecx
	movl	$1717986919, %edx
	movl	%ecx, %eax
	imull	%edx
	sarl	$3, %edx
	movl	%ecx, %eax
	sarl	$31, %eax
	subl	%eax, %edx
	movl	%edx, %eax
	sall	$2, %eax
	addl	%edx, %eax
	sall	$2, %eax
	subl	%eax, %ecx
	movl	%ecx, %edx
	pxor	%xmm0, %xmm0
	cvtsi2sd	%edx, %xmm0
	movsd	LC0(%rip), %xmm1
	divsd	%xmm1, %xmm0
	movd	%xmm0, %rax
	movd	%rax, %xmm0
	call	_exp
	movd	%xmm0, %r13
	movl	$1759218605, %edx
	movl	%ebx, %eax
	imull	%edx
	sarl	$13, %edx
	movl	%ebx, %eax
	sarl	$31, %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	imull	$20000, %ecx, %eax
	movl	%ebx, %ecx
	subl	%eax, %ecx
	movl	$1098962147, %edx
	movl	%ecx, %eax
	imull	%edx
	sarl	$9, %edx
	movl	%ecx, %eax
	sarl	$31, %eax
	subl	%eax, %edx
	movl	%edx, %eax
	pxor	%xmm2, %xmm2
	cvtsi2sd	%eax, %xmm2
	movd	%xmm2, %rax
	movd	%rax, %xmm0
	call	_sin
	movd	%xmm0, %rax
	movd	%r13, %xmm0
	movd	%rax, %xmm3
	addsd	%xmm3, %xmm0
	leaq	0(,%r14,8), %rdx
	leaq	_a.3384(%rip), %rax
	movsd	%xmm0, (%rdx,%rax)
	addl	$1, %ebx
	subl	$1, %r12d
	jmp	L3
L2:
	movsd	LC1(%rip), %xmm0
	movsd	%xmm0, -64(%rbp)
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	movl	$0, %eax
	call	__gfortran_cpu_time_8
	movl	$1, -36(%rbp)
	cmpl	$20000000, -36(%rbp)
	jg	L4
L7:
	movl	-36(%rbp), %eax
	cltq
	subq	$1, %rax
	leaq	0(,%rax,8), %rdx
	leaq	_a.3384(%rip), %rax
	movsd	(%rdx,%rax), %xmm0
	movsd	-64(%rbp), %xmm1
	ucomisd	%xmm1, %xmm0
	jb	L5
	movl	-36(%rbp), %eax
	cltq
	subq	$1, %rax
	leaq	0(,%rax,8), %rdx
	leaq	_a.3384(%rip), %rax
	movsd	(%rdx,%rax), %xmm0
	movsd	%xmm0, -64(%rbp)
L5:
	cmpl	$20000000, -36(%rbp)
	sete	%al
	movzbl	%al, %eax
	addl	$1, -36(%rbp)
	testl	%eax, %eax
	jne	L4
	jmp	L7
L4:
	leaq	-56(%rbp), %rax
	movq	%rax, %rdi
	movl	$0, %eax
	call	__gfortran_cpu_time_8
	leaq	LC2(%rip), %rax
	movq	%rax, -552(%rbp)
	movl	$20, -544(%rbp)
	movl	$128, -560(%rbp)
	movl	$6, -556(%rbp)
	leaq	-560(%rbp), %rax
	movq	%rax, %rdi
	call	__gfortran_st_write
	movsd	-56(%rbp), %xmm0
	movsd	-48(%rbp), %xmm1
	subsd	%xmm1, %xmm0
	movsd	%xmm0, -72(%rbp)
	leaq	-72(%rbp), %rcx
	leaq	-560(%rbp), %rax
	movl	$8, %edx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	__gfortran_transfer_real_write
	leaq	-64(%rbp), %rcx
	leaq	-560(%rbp), %rax
	movl	$8, %edx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	__gfortran_transfer_real_write
	leaq	-560(%rbp), %rax
	movq	%rax, %rdi
	call	__gfortran_st_write_done
	nop
	addq	$528, %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%rbp
LCFI3:
	ret
LFE0:
	.globl _main
_main:
LFB1:
	pushq	%rbp
LCFI4:
	movq	%rsp, %rbp
LCFI5:
	subq	$16, %rsp
	movl	%edi, -4(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-16(%rbp), %rdx
	movl	-4(%rbp), %eax
	movq	%rdx, %rsi
	movl	%eax, %edi
	call	__gfortran_set_args
	leaq	_options.3.3407(%rip), %rsi
	movl	$9, %edi
	call	__gfortran_set_options
	call	_MAIN__
	movl	$0, %eax
	leave
LCFI6:
	ret
LFE1:
	.zerofill __DATA,__bss5,_a.3384,160000000,5
	.const
	.align 5
_options.3.3407:
	.long	68
	.long	1023
	.long	0
	.long	0
	.long	1
	.long	1
	.long	0
	.long	0
	.long	31
	.literal8
	.align 3
LC0:
	.long	0
	.long	1080623104
	.align 3
LC1:
	.long	4294967295
	.long	-1048577
	.section __TEXT,__eh_frame,coalesced,no_toc+strip_static_syms+live_support
EH_frame1:
	.set L$set$0,LECIE1-LSCIE1
	.long L$set$0
LSCIE1:
	.long	0
	.byte	0x1
	.ascii "zR\0"
	.byte	0x1
	.byte	0x78
	.byte	0x10
	.byte	0x1
	.byte	0x10
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.byte	0x90
	.byte	0x1
	.align 3
LECIE1:
LSFDE1:
	.set L$set$1,LEFDE1-LASFDE1
	.long L$set$1
LASFDE1:
	.long	LASFDE1-EH_frame1
	.quad	LFB0-.
	.set L$set$2,LFE0-LFB0
	.quad L$set$2
	.byte	0
	.byte	0x4
	.set L$set$3,LCFI0-LFB0
	.long L$set$3
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$4,LCFI1-LCFI0
	.long L$set$4
	.byte	0xd
	.byte	0x6
	.byte	0x4
	.set L$set$5,LCFI2-LCFI1
	.long L$set$5
	.byte	0x8e
	.byte	0x3
	.byte	0x8d
	.byte	0x4
	.byte	0x8c
	.byte	0x5
	.byte	0x83
	.byte	0x6
	.byte	0x4
	.set L$set$6,LCFI3-LCFI2
	.long L$set$6
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.align 3
LEFDE1:
LSFDE3:
	.set L$set$7,LEFDE3-LASFDE3
	.long L$set$7
LASFDE3:
	.long	LASFDE3-EH_frame1
	.quad	LFB1-.
	.set L$set$8,LFE1-LFB1
	.quad L$set$8
	.byte	0
	.byte	0x4
	.set L$set$9,LCFI4-LFB1
	.long L$set$9
	.byte	0xe
	.byte	0x10
	.byte	0x86
	.byte	0x2
	.byte	0x4
	.set L$set$10,LCFI5-LCFI4
	.long L$set$10
	.byte	0xd
	.byte	0x6
	.byte	0x4
	.set L$set$11,LCFI6-LCFI5
	.long L$set$11
	.byte	0xc
	.byte	0x7
	.byte	0x8
	.align 3
LEFDE3:
	.subsections_via_symbols
