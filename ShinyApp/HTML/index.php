<!DOCTYPE html>
<html lang="en">
	<head>
		<meta charset="UTF-8" />
		<meta name="viewport" content="width=device-width, initial-scale=1.0"> 
		<title>Network Visualization</title>
		<meta name="description" content="Network Visualization Panel" />
		<meta name="keywords" content="Admin, Dashboard, Bootstrap3, Sass, transform, CSS3, HTML5, Web design, UI Design, Responsive Dashboard, Responsive Admin, Admin Theme, Best Admin UI, Bootstrap Theme, Bootstrap, Light weight Admin Dashboard,Light weight, Light weight Admin, Light weight Dashboard" />
		<meta name="author" content="Bootstrap Gallery" />
		<link rel="shortcut icon" href="img/favicon.ico">
	
		<!-- Bootstrap CSS -->
		<link href="css/bootstrap.min.css" rel="stylesheet" media="screen">

		<!-- Main CSS -->
		<link href="css/main.css" rel="stylesheet" media="screen">

		<!-- Font Awesome -->
		<link href="fonts/font-awesome.min.css" rel="stylesheet">

		<!-- Metrize Fonts -->
		<link href="fonts/metrize.css" rel="stylesheet">

		<!-- HTML5 shiv and Respond.js IE8 support of HTML5 elements and media queries -->
		<!--[if lt IE 9]>
			<script src="js/html5shiv.js"></script>
			<script src="js/respond.min.js"></script>
		<![endif]-->

		<style type="text/css">
			iframe {
				width: 100%;
				height: 80vh;
			}

			#loading {
				display: block;
				position: absolute;
				top: 0;
				left: 0;
				z-index: 100;
				width: 100vw;
				height: 150vh;
				background-color: rgba(192, 192, 192, 0.5);
				background-image: url("https://i.stack.imgur.com/MnyxU.gif");
				background-repeat: no-repeat;
				background-position: center;
			}
		</style>

	</head>  

	<body>

		<div id="loading"></div>

		<div id="page">

			<!-- Left sidebar start -->
			<aside id="sidebar">

				<!-- Logo starts -->
				<a href="#" class="logo">
					<!--<img src="img/logo.png" alt="logo">-->
				</a>
				<!-- Logo ends -->

				<!-- Menu start -->
				<div id='menu'>
					<?php

						/**
						 * Define the sort function
						 */
						function sortByIndex($a, $b)
						{
							$a = $a['index'];
							$b = $b['index'];

							if ($a == $b) return 0;
							return ($a < $b) ? -1 : 1;
						}

						// Initialize the objects
						$menuList = array();
						$menuHTML = '';
						$idx = 1;

						// Run through the files in folder and fill the array
						if ($handle = opendir('./networks/')) {
							while (false !== ($file = readdir($handle)))
							{
								if ($file != "." && $file != ".." && strtolower(substr($file, strrpos($file, '.') + 1)) == 'html')
								{
									$menuList[$idx]['pathwayCode'] = substr($file, -10, 5);
									$menuList[$idx]['index'] = intval(explode('_', $file, -1)[0]);
									$idx += 1;
								}
							}
						}

						// Check if the menuList is empty
						if(!empty($menuList)) {
							// Sort the array by the pathway name
							usort($menuList, 'sortByIndex');
						}

						// Print the menus
						foreach ($menuList as $item => $value) { 
							$$menuHTML = '';
							$$menuHTML .= '<ul style="cursor: pointer;">';

							if ($value['index'] == 1) {
								$$menuHTML .= '<li id="menu' . $value['index'] . '" name="menuList" class="highlight">';
							} else {
								$$menuHTML .= '<li id="menu' . $value['index'] . '" name="menuList" >';
							}

							$$menuHTML .= '<a onclick="setNetwork(\'' . $value['pathwayCode'] . '\',\'' . $value['index'] . '\');' .
											'setNetworkAttributes(\'menu' . $value['index'] . '\',\'' . $value['index'] . '\')">';

							$$menuHTML .= '<div class="fs1" aria-hidden="true" data-icon="&#xe0f8;"></div>';
							$$menuHTML .= '<span>' . $value['pathwayCode']  . '</span>';

							$$menuHTML .= '</li>';
							$$menuHTML .= '</ul>';

							echo $$menuHTML;
						}
					?>
				</div>
				<!-- Menu End -->

				<!-- Extras starts -->
				<div class="extras">
					<div class="ex-wrapper">
						<div class="alert alert-info">
							<strong>Important!</strong> This project still in progress, so the visualizations can be modified without warning.
						</div>
					</div>
				</div>
				<!-- Extras ends -->
				
			</aside>
			<!-- Left sidebar end -->

			<!-- Dashboard Wrapper Start -->
			<div class="dashboard-wrapper">

				<!-- Header start -->
				<header>
					<ul class="pull-left" id="left-nav">
						<li class="hidden-lg hidden-md hidden-sm">
							<div class="logo-mob clearfix">
								<h2><div class="fs1" aria-hidden="true" data-icon="&#xe0db;"></div><span>Network Visualization</span></h2>
							</div>
						</li>
						<li>
							<div class="custom-search hidden-sm hidden-xs pull-left">
								<input type="text" class="search-query" placeholder="Search here">
								<i class="fa fa-search"></i>
							</div>
						</li>
					</ul>
					<div class="pull-right"></div>
				</header>
				<!-- Header ends -->

				<!-- Main Container Start -->
				<div class="main-container">

					<!-- Top Bar Starts -->
					<div class="top-bar clearfix">
						<div class="page-title">
							<h4>
								<div class="fs1" aria-hidden="true" data-icon="&#xe0f8;"></div>
								<span id='graphTitle'>Graph title</span>
							</h4>
						</div>
					</div>
					<!-- Top Bar Ends -->

					<!-- Container fluid Starts -->
					<div class="container-fluid">

						<!-- Spacer starts -->
						<div class="spacer-xs">

							<!-- Row start -->
							<div class="row no-gutter">
								<div class="col-md-12 col-sm-12 col-sx-12">
									<div class="panel">
										<div class="panel-heading">
											<h4>&nbsp;</h4>
										</div>
										<div class="panel-body">
											<!-- Netwrok iframe -->
											<iframe id="networkIframe" src="http://getprismatic.com/" frameBorder="0"></iframe>
										</div>
									</div>
								</div>
							</div>
							<!-- Row end -->

						</div>
						<!-- Spacer ends -->

					</div>
					<!-- Container fluid ends -->

				</div>
				<!-- Main Container Start -->

				<!-- Footer Start -->
				<footer>
					Developed by <a href="https://igorabrandao.com.br" target="_blank"> Igor Brandão </a>
				</footer>
				<!-- Footer end -->
				
			</div>
			<!-- Dashboard Wrapper End -->

		</div>

		<!-- jQuery (necessary for Bootstrap's JavaScript plugins) -->
		<script src="js/jquery.js"></script>

		<!-- Custom .js code -->
		<script src="js/custom.js"></script>

		<!-- Include all compiled plugins (below), or include individual files as needed -->
		<script src="js/bootstrap.min.js"></script>

		<!-- Custom JS -->
		<script src="js/custom.js"></script>
	</body>
</html>